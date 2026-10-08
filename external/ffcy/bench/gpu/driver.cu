// Copyright (C) 2026 Jack S. Hale
//
// This file is part of FFCy (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    MIT

// Checks and times an action compiled by MLIR, for any form.
//
// The kernel is a shared library. A case written by bench/cases.py gives the arguments
// and expected result on a few cells, and says which of two kinds of action it is:
//
// - On E-vectors, _mlir_ciface_action(input, input, output), whose memrefs have a
//   leading cell axis. The case's cells are repeated to the problem size.
// - On assembled vectors, _mlir_ciface_action(geometry, x, dofmap, values), whose
//   values on cells the driver sums into y. The case is a small box, and
//   timing runs on a larger box numbered the same way.
//
// With --cells-per-block C, arrays on cells are split into blocks of C cells. Each dof
// sums its entries through the transposed dofmap in CSR form, which is deterministic.
// With --scatter fused the action itself adds its values into y atomically,
// _mlir_ciface_action(geometry, x, dofmap, y).

#include "../case.h"

#include <cuda_runtime.h>
#include <dlfcn.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace
{
using namespace ffcy;

void check(cudaError_t error, const char* what)
{
  if (error != cudaSuccess)
    fail(std::string(what) + ": " + cudaGetErrorString(error));
}

template <class T>
struct DeviceArray
{
  T* data = nullptr;
  std::int64_t size;

  explicit DeviceArray(const std::vector<T>& host) : size(host.size())
  {
    check(cudaMalloc(&data, size * sizeof(T)), "cudaMalloc");
    check(cudaMemcpy(data, host.data(), size * sizeof(T), cudaMemcpyHostToDevice),
          "cudaMemcpy");
  }
  DeviceArray(const DeviceArray&) = delete;
  ~DeviceArray() { cudaFree(data); }

  std::vector<T> host() const
  {
    std::vector<T> h(size);
    check(cudaMemcpy(h.data(), data, size * sizeof(T), cudaMemcpyDeviceToHost),
          "cudaMemcpy");
    return h;
  }

  std::int64_t bytes() const { return size * sizeof(T); }
};

// Milliseconds per call of each function in `steps`, run in turn for at least
// `seconds` in total.
std::vector<double> time(const std::vector<std::function<void()>>& steps, double seconds)
{
  std::vector<cudaEvent_t> events(steps.size() + 1);
  for (cudaEvent_t& e : events)
    check(cudaEventCreate(&e), "cudaEventCreate");
  std::vector<double> ms(steps.size(), 0);
  int calls = 0;
  for (double total = 0; calls < 5 or total < 1000 * seconds; ++calls)
  {
    check(cudaEventRecord(events[0]), "cudaEventRecord");
    for (std::size_t s = 0; s < steps.size(); ++s)
    {
      steps[s]();
      check(cudaEventRecord(events[s + 1]), "cudaEventRecord");
    }
    check(cudaEventSynchronize(events.back()), "cudaEventSynchronize");
    for (std::size_t s = 0; s < steps.size(); ++s)
    {
      float elapsed;
      check(cudaEventElapsedTime(&elapsed, events[s], events[s + 1]), "cudaEventElapsedTime");
      ms[s] += elapsed;
      total += elapsed;
    }
  }
  for (double& m : ms)
    m /= calls;
  return ms;
}

void run_on_cells(void* library, const std::string& prefix, std::int64_t cells_per_block,
                  double num_dofs, double seconds)
{
  auto action = reinterpret_cast<void (*)(void*, void*, void*)>(
      dlsym(library, "_mlir_ciface_action"));
  if (!action)
    fail(dlerror());
  std::ifstream data;
  auto lines = read_case(prefix, data);
  const std::int64_t case_cells = lines["cells"][0];
  std::vector<Shape> shapes;
  for (auto& [word, numbers] : lines)
    if (word.starts_with("input") or word.starts_with("output"))
      shapes.push_back(numbers);
  if (shapes.size() != 3)
    fail("expected an action with two inputs and one output");

  const std::int64_t block = std::max<std::int64_t>(cells_per_block, 1);
  const std::int64_t dofs_per_cell = product(shapes[2]);
  const std::int64_t num_cells
      = (std::llround(num_dofs / dofs_per_cell) + block - 1) / block * block;
  std::vector<std::vector<double>> base;
  for (const Shape& shape : shapes)
    base.push_back(read_array<double>(data, case_cells * product(shape)));
  const std::vector<double> expected = base[2];
  // The output starts as NaN, so a kernel that does not run fails the check.
  std::fill(base[2].begin(), base[2].end(), std::nan(""));

  std::deque<DeviceArray<double>> arrays;
  std::vector<std::vector<std::int64_t>> descriptors;
  for (std::size_t i = 0; i < 3; ++i)
  {
    auto& a = arrays.emplace_back(repeat(base[i], num_cells * product(shapes[i])));
    descriptors.push_back(descriptor(a.data, on_cells(num_cells, cells_per_block, shapes[i])));
  }
  auto apply = [&]
  { action(descriptors[0].data(), descriptors[1].data(), descriptors[2].data()); };

  apply();
  check(cudaDeviceSynchronize(), "action");
  const double error = relative_error(arrays[2].host(), expected);
  double bytes = 0;
  for (const auto& a : arrays)
    bytes += a.bytes();
  report(error, double(num_cells) * dofs_per_cell, time({apply}, seconds), bytes);
}

constexpr int threads_per_block = 256;

__global__ void scatter_csr(const double* values, const std::int32_t* offsets,
                            const std::int32_t* entries, std::int64_t num_dofs, double* y)
{
  const std::int64_t dof = blockIdx.x * std::int64_t{blockDim.x} + threadIdx.x;
  if (dof >= num_dofs)
    return;
  double sum = 0;
  for (std::int32_t k = offsets[dof]; k < offsets[dof + 1]; ++k)
    sum += values[entries[k]];
  y[dof] = sum;
}

std::int64_t num_blocks(std::int64_t size)
{
  return (size + threads_per_block - 1) / threads_per_block;
}

void run_assembled(void* library, const std::string& prefix, std::int64_t cells_per_block,
                   double target_dofs, double seconds, bool fused)
{
  auto action = reinterpret_cast<void (*)(void*, void*, void*, void*)>(
      dlsym(library, "_mlir_ciface_action"));
  if (!action)
    fail(dlerror());
  std::ifstream data;
  auto lines = read_case(prefix, data);
  const Shape geometry_shape = lines["geometry"], values_shape = lines["values"];
  const int d = values_shape.size();
  const std::int64_t degree = values_shape[0] - 1, ndofs = product(values_shape);
  const std::int64_t width = std::int64_t{1} << d;
  const std::vector<std::int64_t>& lattice = lines["lattice"];

  auto apply = [&](std::int64_t num_cells, std::int64_t num_dofs,
                   const std::vector<double>& geometry, const std::vector<double>& x,
                   const std::vector<std::int32_t>& dofmap,
                   const std::vector<std::int32_t>& transpose, double seconds)
  {
    DeviceArray g(geometry), v(x);
    DeviceArray m(dofmap);
    const auto [host_offsets, host_entries] = compress(transpose, width);
    DeviceArray offsets(host_offsets), entries(host_entries);
    DeviceArray values(std::vector<double>(num_cells * ndofs, std::nan("")));
    DeviceArray y(std::vector<double>(num_dofs, std::nan("")));
    auto gd = descriptor(g.data, on_cells(num_cells, cells_per_block, geometry_shape));
    auto vd = descriptor(v.data, {num_dofs});
    auto md = descriptor(m.data, on_cells(num_cells, cells_per_block, {ndofs}));
    auto cd = descriptor(values.data, on_cells(num_cells, cells_per_block, values_shape));
    auto yd = descriptor(y.data, {num_dofs});
    std::vector<std::function<void()>> steps{
        [&] { action(gd.data(), vd.data(), md.data(), cd.data()); },
        [&]
        {
          scatter_csr<<<num_blocks(num_dofs), threads_per_block>>>(
              values.data, offsets.data, entries.data, num_dofs, y.data);
        }};
    if (fused)
      steps = {[&] { check(cudaMemsetAsync(y.data, 0, y.bytes()), "cudaMemsetAsync"); },
               [&] { action(gd.data(), vd.data(), md.data(), yd.data()); }};
    for (auto& step : steps)
      step();
    check(cudaDeviceSynchronize(), "action");
    std::vector<double> ms;
    if (seconds > 0)
      ms = time(steps, seconds);
    // Atomic adds zero y then update it. The scatter also writes and reads the values,
    // and reads the transpose.
    const double scatter_bytes
        = fused ? 3 * y.bytes()
                : 2 * values.bytes() + offsets.bytes() + entries.bytes() + y.bytes();
    const double bytes = g.bytes() + v.bytes() + m.bytes() + scatter_bytes;
    return std::tuple{y.host(), ms, bytes};
  };

  // Check on the case's box.
  const std::int64_t case_cells = lines["cells"][0], case_dofs = lines["dofs"][0];
  const auto geometry = read_array<double>(data, case_cells * product(geometry_shape));
  const auto x = read_array<double>(data, case_dofs);
  const auto dofmap = read_array<std::int32_t>(data, case_cells * ndofs);
  const auto transpose = read_array<std::int32_t>(data, case_dofs * width);
  const auto expected = read_array<double>(data, case_dofs);
  const auto [y, unused, unused_bytes]
      = apply(case_cells, case_dofs, geometry, x, dofmap, transpose, 0);
  const double error = relative_error(y, expected);

  // Time on a box with about target_dofs dofs, whose cells split into blocks of 64.
  const std::int64_t step = d == 3 ? 4 : 8;
  std::int64_t n = std::llround((std::pow(target_dofs, 1.0 / d) - 1) / degree);
  n = std::max(step, (n + step - 1) / step * step);
  const Box box(n, degree, d, lattice);
  const auto [unused_y, ms, bytes] = apply(
      box.num_cells, box.num_dofs, repeat(geometry, box.num_cells * product(geometry_shape)),
      repeat(x, box.num_dofs), box.dofmap, box.transpose, seconds);
  report(error, box.num_dofs, ms, bytes,
         fused ? std::vector<std::string>{"zero", "action"}
               : std::vector<std::string>{"action", "scatter"});
}
} // namespace

int main(int argc, char** argv)
{
  if (argc < 3)
    fail("usage: ffcy_driver KERNEL.so CASE_PREFIX [--cells-per-block C] [--num-dofs N] "
         "[--seconds S] [--scatter fused]");
  std::int64_t cells_per_block = 0;
  bool fused = false;
  double num_dofs = 4e6;
  double seconds = 3;
  for (int i = 3; i + 1 < argc; i += 2)
  {
    if (!std::strcmp(argv[i], "--cells-per-block"))
      cells_per_block = std::atoll(argv[i + 1]);
    else if (!std::strcmp(argv[i], "--num-dofs"))
      num_dofs = std::atof(argv[i + 1]);
    else if (!std::strcmp(argv[i], "--seconds"))
      seconds = std::atof(argv[i + 1]);
    else if (!std::strcmp(argv[i], "--scatter") and !std::strcmp(argv[i + 1], "fused"))
      fused = true;
    else
      fail(std::string("unknown option ") + argv[i]);
  }

  void* library = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
  if (!library)
    fail(dlerror());
  std::ifstream data;
  if (read_case(argv[2], data).count("dofs"))
    run_assembled(library, argv[2], cells_per_block, num_dofs, seconds, fused);
  else
    run_on_cells(library, argv[2], cells_per_block, num_dofs, seconds);
}
