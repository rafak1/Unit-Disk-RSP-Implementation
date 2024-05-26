cmake_minimum_required(VERSION 3.0)
project(alg_2)

add_library(alg_2 INTERFACE)

target_include_directories(alg_2 INTERFACE include)

target_link_libraries(alg_2 INTERFACE grid)
target_link_libraries(alg_2 INTERFACE cs)