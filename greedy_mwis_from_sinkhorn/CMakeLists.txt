cmake_minimum_required(VERSION 3.14)
project(CombinatorialSolver CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(BSA_GREEDY main.cpp)

# Add include path where you put json.hpp
target_include_directories(BSA_GREEDY PRIVATE ${CMAKE_SOURCE_DIR}/external/)
