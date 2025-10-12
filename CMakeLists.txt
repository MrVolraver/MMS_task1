cmake_minimum_required(VERSION 3.31)
project(mms1)

set(CMAKE_CXX_STANDARD 20)

add_executable(main main.cpp
)

add_library(RungeKutta SHARED RungeKutta.cpp)

add_library(Adams_Bashforth_Moulton SHARED Adams_Bashforth_Moulton.cpp)

add_library(Milne_Simpson SHARED Milne_Simpson.cpp)