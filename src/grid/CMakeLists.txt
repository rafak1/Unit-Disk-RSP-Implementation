cmake_minimum_required(VERSION 3.0)
project(grid)

add_library(grid INTERFACE)

target_include_directories(grid INTERFACE include)

export(TARGETS grid FILE gridConfig.cmake)