# We don't do anything fancy here
cmake_minimum_required(VERSION 2.6)
set(CMAKE_BUILD_TYPE Debug)

# Give your project a name
PROJECT (baofisher)

# Includes 
include_directories("$ENV{NPTOOLS_DIR}/include")
link_directories("$ENV{NPTOOLS_DIR}/lib")


add_library(baofisher SHARED bao_forecast.c)

add_executable(unittest unittest.cpp)
target_link_libraries(unittest gtest gtest_main baofisher)



