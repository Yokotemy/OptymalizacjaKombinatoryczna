cmake_minimum_required(VERSION 4.0)
project(LocalSearch)

set(CMAKE_CXX_STANDARD 20)

add_executable(LocalSearch main.cpp
        TSP.cpp
        TSP.h)
