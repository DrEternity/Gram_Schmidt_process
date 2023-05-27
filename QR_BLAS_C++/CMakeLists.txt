cmake_minimum_required(VERSION 3.24)
project(QR)

enable_language(Fortran)

set(GCC_COVERAGE_COMPILE_FLAGS "-flto -Ofast -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} -lcblas")

find_package(BLAS REQUIRED)

add_executable(QR QR.cpp)

