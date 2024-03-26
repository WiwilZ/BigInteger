cmake_minimum_required(VERSION 3.22.1)
project(Integer)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -m64")

add_executable(Integer main.cpp
        Utils.h
        BigInteger.h
)
