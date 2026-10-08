cmake_minimum_required(VERSION 3.24)
project(ffcy_cpu_plugin LANGUAGES CXX)

include(../add_plugin.cmake)
add_ffcy_plugin(FFCyCpuPlugin plugin.cpp slice_constant_indices.cpp
                tile_stages.cpp pipeline.cpp)
