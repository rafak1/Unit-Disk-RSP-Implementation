cmake_minimum_required(VERSION 3.0)
project(cs_algorithm)

option(BUILD_SHARED_LIBS ON)

include_directories(grid)
include_directories(cs)

add_subdirectory(cs)
add_subdirectory(grid)
add_executable(main main.cpp)

target_link_libraries( main
    INTERFACE 
        grid
        cs
)