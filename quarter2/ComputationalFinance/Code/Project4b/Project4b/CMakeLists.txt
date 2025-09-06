cmake_minimum_required(VERSION 3.25)
project(Project4b)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(Project4b
        greeks.cpp
        greeks.hpp
        halton_sequence.cpp
        halton_sequence.hpp
        main.cpp
        max.hpp
        option.cpp
        option.hpp)
