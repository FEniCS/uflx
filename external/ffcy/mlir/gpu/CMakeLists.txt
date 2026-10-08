cmake_minimum_required(VERSION 3.24)
project(ffcy_gpu_plugin LANGUAGES CXX)

include(../add_plugin.cmake)
add_ffcy_plugin(FFCyGpuPlugin plugin.cpp workgroup_allocations.cpp
                atomic_scatter.cpp pencil_schedule.cpp pipeline.cpp)
