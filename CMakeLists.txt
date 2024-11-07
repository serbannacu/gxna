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

# Main target

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

# Unit tests

enable_testing()
find_package(GTest)

if(GTest_FOUND)
  add_executable(test_distribution test/test_distribution.cpp cpp/Distribution.cpp)
  add_executable(test_statistics test/test_statistics.cpp cpp/Statistics.cpp)

  target_link_libraries(test_distribution GTest::gtest GTest::gtest_main)
  target_link_libraries(test_statistics GTest::gtest GTest::gtest_main)

  include(GoogleTest)
  gtest_discover_tests(test_distribution)
  gtest_discover_tests(test_statistics)
endif()

# Regression tests

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
  NAME reg_gxna_t
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND test/regression_test.sh 102 -algoType GXNA -depth 15 -draw T
)

add_test(
  NAME reg_gxna_f
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND test/regression_test.sh 103 -algoType GXNA -depth 15 -draw T -test group
)
