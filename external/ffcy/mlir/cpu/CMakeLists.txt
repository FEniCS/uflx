# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

cmake_minimum_required(VERSION 3.24)
project(ffcy_cpu_plugin LANGUAGES CXX)

include(../add_plugin.cmake)
add_ffcy_plugin(FFCyCpuPlugin plugin.cpp slice_constant_indices.cpp
                tile_stages.cpp pipeline.cpp)
