# We don't do anything fancy here
cmake_minimum_required(VERSION 2.6)
set(CMAKE_BUILD_TYPE Debug)

# Give your project a name
PROJECT (gtest_example)

# Includes 
include_directories("$ENV{NPTOOLS_DIR}/include")
link_directories("$ENV{NPTOOLS_DIR}/lib")


add_executable(eg1_test eg1_test.cpp eg1.cpp)
target_link_libraries(eg1_test gtest gtest_main)



