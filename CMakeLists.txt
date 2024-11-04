cmake_minimum_required(VERSION 3.10...3.30)

project(gxna
  VERSION 3.0
  DESCRIPTION "Gene Expression Network Analysis"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_compile_options(-Wall -Wextra)
include_directories(include)

add_executable(gxna
  cpp/gxna.cpp
  cpp/Args.cpp
  cpp/Distribution.cpp
  cpp/Experiment.cpp
  cpp/GeneNetwork.cpp
  cpp/RandomPermutation.cpp
  cpp/Statistics.cpp
)

target_compile_options(gxna PRIVATE -O3)

add_executable(test_distribution
  test/test_distribution.cpp
  cpp/Distribution.cpp
)

add_executable(test_statistics
  test/test_statistics.cpp
  cpp/Statistics.cpp
)

enable_testing()

add_test(
  NAME reg_solo
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND test/regression_test.sh 100 -algoType Basic -radius 0
)

add_test(
  NAME reg_basic
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND test/regression_test.sh 101 -algoType Basic -radius 1
)

add_test(
  NAME reg_gxna
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND test/regression_test.sh 102 -algoType GXNA -depth 15 -draw T
)

add_test(
  NAME unit_dist
  COMMAND test_distribution
)

add_test(
  NAME unit_stat
  COMMAND test_statistics
)
