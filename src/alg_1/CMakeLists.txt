cmake_minimum_required(VERSION 3.0)
project(alg_1)

add_library(alg_1 INTERFACE)

target_include_directories(alg_1 INTERFACE include)

target_link_libraries(alg_1 INTERFACE grid)
target_link_libraries(alg_1 INTERFACE cs)