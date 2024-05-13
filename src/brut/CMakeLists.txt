cmake_minimum_required(VERSION 3.0)
project(brut)

add_library(brut INTERFACE)

target_include_directories(brut INTERFACE include)

export(TARGETS brut FILE brutConfig.cmake)