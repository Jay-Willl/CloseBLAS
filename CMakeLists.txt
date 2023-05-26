cmake_minimum_required(VERSION 3.25)
project(CloseBLAS)

set(CMAKE_CXX_STANDARD 23)

include_directories(/opt/homebrew/Cellar/openblas/0.3.23/include)

link_libraries("/opt/homebrew/Cellar/openblas/0.3.23/lib/libopenblas.dylib")

add_executable(CloseBLAS_Test test.cpp)


