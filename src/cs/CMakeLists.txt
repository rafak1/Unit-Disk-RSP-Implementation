cmake_minimum_required(VERSION 3.0)
project(cs)

add_library(cs INTERFACE)

target_include_directories(cs INTERFACE include)

target_link_libraries(cs INTERFACE grid)

export(TARGETS cs FILE gridConfig.cmake)