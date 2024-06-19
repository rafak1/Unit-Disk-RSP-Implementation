cmake_minimum_required(VERSION 3.0)
project(cs_algorithm)

option(BUILD_SHARED_LIBS ON)

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

# Link the executable to CGAL and third-party libraries

include_directories(grid)
include_directories(cs)
include_directories(alg_1)
include_directories(brut)
include_directories(generator)

add_subdirectory(cs)
add_subdirectory(grid)
add_subdirectory(alg_1)
add_subdirectory(brut)
add_subdirectory(generator)


target_link_libraries(main PRIVATE CGAL::CGAL )