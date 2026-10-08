// Copyright (C) 2026 Jack S. Hale
//
// This file is part of FFCy (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    MIT

// Reads, repeats and checks the cases written by bench/emit.py, for the drivers of
// every target.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace ffcy
{
using Shape = std::vector<std::int64_t>;

[[noreturn]] inline void fail(const std::string& message)
{
  std::fprintf(stderr, "%s\n", message.c_str());
  std::exit(1);
}

inline std::int64_t product(const Shape& shape)
{
  return std::accumulate(shape.begin(), shape.end(), std::int64_t{1}, std::multiplies{});
}

// The layout of StridedMemRefType<T, rank>, for a rank known at run time:
// allocated and aligned pointers, offset, sizes, then strides of a row-major array.
inline std::vector<std::int64_t> descriptor(const void* data, const Shape& shape)
{
  const auto pointer = reinterpret_cast<std::int64_t>(data);
  std::vector<std::int64_t> d{pointer, pointer, 0};
  d.insert(d.end(), shape.begin(), shape.end());
  Shape strides(shape.size());
  std::int64_t stride = 1;
  for (std::size_t i = shape.size(); i-- > 0;)
  {
    strides[i] = stride;
    stride *= shape[i];
  }
  d.insert(d.end(), strides.begin(), strides.end());
  return d;
}

// The shape of an array on cells, split into blocks if cells_per_block > 0.
inline Shape on_cells(std::int64_t num_cells, std::int64_t cells_per_block,
                      const Shape& trailing)
{
  Shape shape = cells_per_block > 0 ? Shape{num_cells / cells_per_block, cells_per_block}
                                     : Shape{num_cells};
  shape.insert(shape.end(), trailing.begin(), trailing.end());
  return shape;
}

// `base` repeated to `size` entries.
template <class T>
std::vector<T> repeat(const std::vector<T>& base, std::int64_t size)
{
  std::vector<T> out(size);
  for (std::int64_t i = 0; i < size; ++i)
    out[i] = base[i % base.size()];
  return out;
}

// The largest difference from the expected values, relative to their largest value.
// NaN means the kernel did not write its output.
inline double relative_error(const std::vector<double>& y, const std::vector<double>& expected)
{
  double error = 0, scale = 0;
  for (std::size_t i = 0; i < y.size(); ++i)
  {
    const double e = expected[i % expected.size()];
    error = std::max(error, std::isnan(y[i]) ? INFINITY : std::abs(y[i] - e));
    scale = std::max(scale, std::abs(e));
  }
  return error / scale;
}

// A case file's lines, by their first word.
inline std::map<std::string, std::vector<std::int64_t>> read_case(const std::string& prefix,
                                                           std::ifstream& data)
{
  std::ifstream text(prefix + ".case");
  if (!text)
    fail("cannot read " + prefix + ".case");
  data.open(prefix + ".bin", std::ios::binary);
  std::map<std::string, std::vector<std::int64_t>> lines;
  std::string line, word;
  for (int i = 0; std::getline(text, line); ++i)
  {
    std::istringstream words(line);
    words >> word;
    // Inputs and outputs repeat, so they keep their order.
    if (word == "input" or word == "output")
      word += std::to_string(i);
    auto& numbers = lines[word];
    for (std::int64_t n; words >> n;)
      numbers.push_back(n);
  }
  return lines;
}

template <class T>
std::vector<T> read_array(std::ifstream& data, std::int64_t size)
{
  std::vector<T> v(size);
  data.read(reinterpret_cast<char*>(v.data()), size * sizeof(T));
  if (!data)
    fail("case data is too short");
  return v;
}

inline void report(double error, double dofs, const std::vector<double>& ms, double bytes,
            const std::vector<std::string>& steps = {})
{
  const double total = std::accumulate(ms.begin(), ms.end(), 0.0);
  std::printf("error %.1e, %.2e dofs, %.3f ms, %.2f GDoF/s, %.0f GB/s", error, dofs, total,
              dofs / total / 1e6, bytes / total / 1e6);
  for (std::size_t s = 0; s < steps.size(); ++s)
    std::printf(", %s %.3f ms", steps[s].c_str(), ms[s]);
  std::printf("\n");
}

// A box of cells numbered as bench/emit.py does: cells in row-major order, and each
// cell's dofs at its lattice positions on the global lattice of degree * n + 1 points
// along each axis.
struct Box
{
  std::int64_t num_cells, num_dofs;
  std::vector<std::int32_t> dofmap, transpose;

  Box(std::int64_t n, std::int64_t degree, int d, const std::vector<std::int64_t>& lattice)
  {
    const std::int64_t ndofs = lattice.size() / d, size = degree * n + 1;
    num_cells = std::llround(std::pow(n, d));
    num_dofs = std::llround(std::pow(size, d));
    dofmap.resize(num_cells * ndofs);
    for (std::int64_t c = 0; c < num_cells; ++c)
      for (std::int64_t i = 0; i < ndofs; ++i)
      {
        std::int64_t dof = 0, rest = c, cell_stride = num_cells, stride = num_dofs;
        for (int a = 0; a < d; ++a)
        {
          cell_stride /= n;
          stride /= size;
          dof += ((rest / cell_stride) * degree + lattice[i * d + a]) * stride;
          rest %= cell_stride;
        }
        dofmap[c * ndofs + i] = dof;
      }
    // A dof is shared by at most 2^d cells, and unused entries are -1.
    const std::int64_t width = std::int64_t{1} << d;
    transpose.assign(num_dofs * width, -1);
    std::vector<std::int64_t> filled(num_dofs, 0);
    for (std::int64_t e = 0; e < num_cells * ndofs; ++e)
      transpose[dofmap[e] * width + filled[dofmap[e]]++] = e;
  }
};
// The ELL transpose without its padding, as offsets into the entries of each dof.
inline std::pair<std::vector<std::int32_t>, std::vector<std::int32_t>>
compress(const std::vector<std::int32_t>& transpose, std::int64_t width)
{
  std::vector<std::int32_t> offsets{0}, entries;
  for (std::size_t dof = 0; dof < transpose.size() / width; ++dof)
  {
    for (std::int64_t k = 0; k < width; ++k)
      if (const std::int32_t e = transpose[dof * width + k]; e >= 0)
        entries.push_back(e);
    offsets.push_back(entries.size());
  }
  return {offsets, entries};
}
} // namespace ffcy
