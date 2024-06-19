# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( praca_roczna )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE "Release")


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package


# Creating entries for target: praca_roczna
# ############################

add_to_cached_list( CGAL_EXECUTABLE_TARGETS praca_roczna )

include_directories(src/grid)
include_directories(src/cs)
include_directories(src/alg_1)
include_directories(src/brut)
include_directories(src/generator)

add_subdirectory(src/cs)
add_subdirectory(src/grid)
add_subdirectory(src/alg_1)
add_subdirectory(src/brut)
add_subdirectory(src/generator)

add_executable(test_generator src/generator.cpp)

target_link_libraries( test_generator
    INTERFACE 
        generator
)

add_executable(alg_1_solver src/alg_1.cpp)

target_link_libraries( alg_1_solver
    INTERFACE 
        alg_1
        grid
)

add_executable(alg_2_solver src/alg_2.cpp)

target_link_libraries( alg_2_solver
    INTERFACE 
        alg_2
        grid
)

target_link_libraries(alg_2_solver PRIVATE CGAL::CGAL )

add_executable(brut_solver src/brut.cpp)

target_link_libraries( brut_solver
    INTERFACE 
        brut
        grid
)


