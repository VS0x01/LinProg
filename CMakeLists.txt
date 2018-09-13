cmake_minimum_required(VERSION 3.12)
project(LinProg)

set(CMAKE_CXX_STANDARD 14)

include_directories(inc)

add_executable(LinProg
        inc/GaussJordan.h
        src/GaussJordan.cpp
        main.cpp)
