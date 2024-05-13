cmake_minimum_required(VERSION 3.0)
project(generator)

add_library(generator INTERFACE)

target_include_directories(generator INTERFACE include)

target_link_libraries(generator INTERFACE grid)