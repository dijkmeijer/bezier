cmake_minimum_required (VERSION 2.8.11)
project (dxflib)

add_library( dxflib STATIC dl_dxf.cpp dl_writer_ascii.cpp )
