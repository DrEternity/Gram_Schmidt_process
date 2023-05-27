cmake_minimum_required(VERSION 3.16)
project(QR)

find_package(LAPACK REQUIRED)


set(GCC_COVERAGE_COMPILE_FLAGS "-flto -Ofast -march=native -llapack")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")


add_executable(QR QR.cpp)

