cmake_minimum_required(VERSION 3.10...3.30)

project(gxna
  VERSION 3.0
  DESCRIPTION "Gene Expression Network Analysis"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(gxna
  cpp/gxna.cpp
  cpp/Args.cpp
  cpp/Distribution.cpp
  cpp/Experiment.cpp
  cpp/GeneNetwork.cpp
  cpp/RandomPermutation.cpp
  cpp/Statistics.cpp
)

target_include_directories(gxna PUBLIC include)
target_compile_options(gxna PRIVATE -Wall -Wextra -g -O3)
