cmake_minimum_required(VERSION 3.24)
project(QR)

enable_language(Fortran)
find_package(LAPACK REQUIRED)


set(GCC_COVERAGE_COMPILE_FLAGS "-flto -Ofast -march=native -llapack")
set(CMAKE_Fortran_FLAGS "${GCC_COVERAGE_COMPILE_FLAGS}")


add_executable(QR QR.f90)

