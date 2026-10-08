// Copyright (C) 2026 Jack S. Hale
//
// This file is part of FFCy (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    MIT

// Checks and times an action compiled by MLIR on CPU cores, for any form.
//
// The kernel is a shared library whose memrefs hold cells interleaved in blocks of C,
// given by --cells-per-block C. An array of shape [cells, ...] is [cells / C, ..., C],
// with each block's cells in the last axis. A case written by bench/cases.py gives the
// arguments and expected result, and says which of two kinds of action it is:
//
// - On E-vectors, _mlir_ciface_action(input, input, output). The case's cells are
//   repeated to the problem size.
// - On assembled vectors, _mlir_ciface_action(geometry, x, dofmap, values), whose
//   values on cells the driver sums into y through the transposed dofmap. The case
//   is a small box, and timing runs on a larger box numbered the same way.
//
// Each OpenMP thread runs the kernel on its own contiguous range of blocks, and
// first touches that range, so its data is on its own NUMA node. Pin threads with
// OMP_PLACES and OMP_PROC_BIND.

#include "../case.h"

#include <dlfcn.h>
#include <omp.h>

#include <chrono>
#include <cstring>
#include <memory>

namespace
{
using namespace ffcy;

// The blocks that thread t of `threads` runs and first touches.
std::pair<std::int64_t, std::int64_t> blocks_of(int t, int threads, std::int64_t blocks)
{
  return {blocks * t / threads, blocks * (t + 1) / threads};
}

template <class T>
using Array = std::unique_ptr<T[], decltype(&std::free)>;

template <class T>
Array<T> allocate(std::int64_t size)
{
  return Array<T>(static_cast<T*>(std::aligned_alloc(64, size * sizeof(T))), std::free);
}

// An array on cells interleaved in blocks of `lanes`, whose cell c holds `size`
// values of cell c % case_cells of `base`. Each thread fills its own blocks.
template <class T>
Array<T> interleave(const std::vector<T>& base, std::int64_t case_cells,
                    std::int64_t num_cells, std::int64_t size, std::int64_t lanes)
{
  Array<T> out = allocate<T>(num_cells * size);
#pragma omp parallel
  {
    auto [first, last] = blocks_of(omp_get_thread_num(), omp_get_num_threads(),
                                    num_cells / lanes);
    for (std::int64_t c = first * lanes; c < last * lanes; ++c)
      for (std::int64_t i = 0; i < size; ++i)
        out[((c / lanes) * size + i) * lanes + c % lanes] = base[(c % case_cells) * size + i];
  }
  return out;
}

// A copy of `values` that threads first touch in the static schedule of a parallel
// loop over its entries, as the scatter's loop over dofs does.
template <class T>
Array<T> spread(const std::vector<T>& values)
{
  Array<T> out = allocate<T>(values.size());
  const std::int64_t size = values.size();
#pragma omp parallel for schedule(static)
  for (std::int64_t i = 0; i < size; ++i)
    out[i] = values[i];
  return out;
}

// The offset of value i of cell c in an array interleaved in blocks of `lanes`.
std::int64_t interleaved(std::int64_t c, std::int64_t i, std::int64_t size, std::int64_t lanes)
{
  return ((c / lanes) * size + i) * lanes + c % lanes;
}

// The values of each cell in turn, from an array on cells interleaved in blocks.
std::vector<double> deinterleave(const double* a, std::int64_t num_cells,
                                 std::int64_t size, std::int64_t lanes)
{
  std::vector<double> out(num_cells * size);
  for (std::int64_t c = 0; c < num_cells; ++c)
    for (std::int64_t i = 0; i < size; ++i)
      out[c * size + i] = a[interleaved(c, i, size, lanes)];
  return out;
}

// Milliseconds per call of each function in `steps`, run in turn for at least
// `seconds` in total.
std::vector<double> time(const std::vector<std::function<void()>>& steps, double seconds)
{
  using clock = std::chrono::steady_clock;
  std::vector<double> ms(steps.size(), 0);
  int calls = 0;
  for (double total = 0; calls < 5 or total < 1000 * seconds; ++calls)
    for (std::size_t s = 0; s < steps.size(); ++s)
    {
      const auto start = clock::now();
      steps[s]();
      const double elapsed
          = std::chrono::duration<double, std::milli>(clock::now() - start).count();
      ms[s] += elapsed;
      total += elapsed;
    }
  for (double& m : ms)
    m /= calls;
  return ms;
}

// Descriptor of the blocks [first, last) of an array on cells with `trailing` shape.
template <class T>
std::vector<std::int64_t> blocks_descriptor(const T* data, std::int64_t first,
                                            std::int64_t last, const Shape& trailing,
                                            std::int64_t lanes)
{
  Shape shape{last - first};
  shape.insert(shape.end(), trailing.begin(), trailing.end());
  shape.push_back(lanes);
  return descriptor(data + first * product(trailing) * lanes, shape);
}

void run_on_cells(void* library, const std::string& prefix, std::int64_t lanes,
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

  const std::int64_t dofs_per_cell = product(shapes[2]);
  const std::int64_t num_cells
      = (std::llround(num_dofs / dofs_per_cell) + lanes - 1) / lanes * lanes;
  const std::int64_t num_blocks = num_cells / lanes;
  std::vector<Array<double>> arrays;
  std::vector<double> expected;
  for (std::size_t i = 0; i < 3; ++i)
  {
    auto base = read_array<double>(data, case_cells * product(shapes[i]));
    // The output starts as NaN, so a kernel that does not run fails the check.
    if (i == 2)
    {
      expected = base;
      std::fill(base.begin(), base.end(), std::nan(""));
    }
    arrays.push_back(interleave(base, case_cells, num_cells, product(shapes[i]), lanes));
  }
  // Each thread runs its blocks through descriptors of its own part of each array.
  auto apply = [&]
  {
#pragma omp parallel
    {
      auto [first, last]
          = blocks_of(omp_get_thread_num(), omp_get_num_threads(), num_blocks);
      if (last > first)
      {
        std::vector<std::vector<std::int64_t>> descriptors;
        for (std::size_t i = 0; i < 3; ++i)
          descriptors.push_back(
              blocks_descriptor(arrays[i].get(), first, last, shapes[i], lanes));
        action(descriptors[0].data(), descriptors[1].data(), descriptors[2].data());
      }
    }
  };

  apply();
  const double error
      = relative_error(deinterleave(arrays[2].get(), num_cells, dofs_per_cell, lanes), expected);
  double bytes = 0;
  for (const Shape& shape : shapes)
    bytes += num_cells * product(shape) * sizeof(double);
  report(error, double(num_cells) * dofs_per_cell, time({apply}, seconds), bytes);
}

void run_assembled(void* library, const std::string& prefix, std::int64_t lanes,
                   double target_dofs, double seconds)
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

  // The action, then the scatter over dofs, on a box of `num_cells` cells whose
  // geometry repeats that of the case's cells.
  auto apply = [&](std::int64_t num_cells, std::int64_t num_dofs,
                   const std::vector<double>& geometry, std::int64_t geometry_cells,
                   const std::vector<double>& x, const std::vector<std::int32_t>& dofmap,
                   const std::vector<std::int32_t>& transpose, double seconds)
  {
    const std::int64_t num_blocks = num_cells / lanes;
    auto g = interleave(geometry, geometry_cells, num_cells, product(geometry_shape), lanes);
    auto m = interleave(dofmap, num_cells, num_cells, ndofs, lanes);
    auto values = allocate<double>(num_cells * ndofs);
    auto y = spread(std::vector<double>(num_dofs, std::nan("")));
    auto xs = spread(x);
    auto vd = descriptor(xs.get(), {num_dofs});
    // The transpose holds values in cell order, which the CSR form keeps interleaved.
    auto [host_offsets, host_entries] = compress(transpose, width);
    for (std::int32_t& e : host_entries)
      e = interleaved(e / ndofs, e % ndofs, ndofs, lanes);
    auto offsets = spread(host_offsets);
    auto entries = spread(host_entries);

    std::vector<std::function<void()>> steps{
        [&]
        {
#pragma omp parallel
          {
            auto [first, last]
                = blocks_of(omp_get_thread_num(), omp_get_num_threads(), num_blocks);
            if (last > first)
            {
              auto gd = blocks_descriptor(g.get(), first, last, geometry_shape, lanes);
              auto md = blocks_descriptor(m.get(), first, last, {ndofs}, lanes);
              auto cd = blocks_descriptor(values.get(), first, last, values_shape, lanes);
              action(gd.data(), vd.data(), md.data(), cd.data());
            }
          }
        },
        [&]
        {
#pragma omp parallel for schedule(static)
          for (std::int64_t dof = 0; dof < num_dofs; ++dof)
          {
            double sum = 0;
            for (std::int32_t k = offsets[dof]; k < offsets[dof + 1]; ++k)
              sum += values[entries[k]];
            y[dof] = sum;
          }
        }};
    for (auto& step : steps)
      step();
    std::vector<double> ms;
    if (seconds > 0)
      ms = time(steps, seconds);
    const double bytes
        = sizeof(double)
              * (num_cells * (product(geometry_shape) + 2 * ndofs) + 2 * num_dofs)
          + sizeof(std::int32_t)
                * (num_cells * ndofs + host_offsets.size() + host_entries.size());
    return std::tuple{std::vector<double>(y.get(), y.get() + num_dofs), ms, bytes};
  };

  // Check on the case's box.
  const std::int64_t case_cells = lines["cells"][0], case_dofs = lines["dofs"][0];
  const auto geometry = read_array<double>(data, case_cells * product(geometry_shape));
  const auto x = read_array<double>(data, case_dofs);
  const auto dofmap = read_array<std::int32_t>(data, case_cells * ndofs);
  const auto transpose = read_array<std::int32_t>(data, case_dofs * width);
  const auto expected = read_array<double>(data, case_dofs);
  if (case_cells % lanes)
    fail("the case's cells do not split into blocks of " + std::to_string(lanes));
  const auto [y, unused, unused_bytes]
      = apply(case_cells, case_dofs, geometry, case_cells, x, dofmap, transpose, 0);
  const double error = relative_error(y, expected);

  // Time on a box with about target_dofs dofs, whose cells split into blocks.
  const std::int64_t step = 4;
  std::int64_t n = std::llround((std::pow(target_dofs, 1.0 / d) - 1) / degree);
  n = std::max(step, (n + step - 1) / step * step);
  const Box box(n, degree, d, lattice);
  if (box.num_cells % lanes)
    fail("the box's cells do not split into blocks of " + std::to_string(lanes));
  const auto [unused_y, ms, bytes]
      = apply(box.num_cells, box.num_dofs, geometry, case_cells, repeat(x, box.num_dofs),
              box.dofmap, box.transpose, seconds);
  report(error, box.num_dofs, ms, bytes, {"action", "scatter"});
}
} // namespace

int main(int argc, char** argv)
{
  if (argc < 3)
    fail("usage: ffcy_driver KERNEL.so CASE_PREFIX [--cells-per-block C] [--num-dofs N] "
         "[--seconds S]");
  std::int64_t cells_per_block = 0;
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
    else
      fail(std::string("unknown option ") + argv[i]);
  }

  if (cells_per_block < 1)
    fail("CPU kernels interleave cells, so --cells-per-block must be given");
  void* library = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
  if (!library)
    fail(dlerror());
  std::ifstream data;
  if (read_case(argv[2], data).count("dofs"))
    run_assembled(library, argv[2], cells_per_block, num_dofs, seconds);
  else
    run_on_cells(library, argv[2], cells_per_block, num_dofs, seconds);
}
