cmake_minimum_required (VERSION 2.8.11)

project (dxf_bezier)

include_directories( ../src)
include_directories( ../bezier)

add_library(bezier STATIC bezier.cpp)
add_executable( dxftest main.cpp)
target_link_libraries( dxftest dxflib bezier)
