cmake_minimum_required(VERSION 3.0)
project(sort)

add_library(sort INTERFACE)

target_include_directories(sort INTERFACE include)

target_link_libraries(sort INTERFACE grid)
target_link_libraries(sort INTERFACE cs)