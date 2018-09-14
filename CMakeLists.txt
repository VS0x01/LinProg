cmake_minimum_required(VERSION 3.12)
project(LinProg)

set(CMAKE_CXX_STANDARD 14)

include_directories(inc)

add_executable(LinProg
        inc/GaussJordan.h
        inc/Simplex.h
        src/GaussJordan.cpp
        src/Simplex.cpp
        main.cpp)
