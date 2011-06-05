# We don't do anything fancy here
cmake_minimum_required(VERSION 2.6)
set(CMAKE_BUILD_TYPE Debug)

# Give your project a name
PROJECT (HewettWild)

# Includes 
include_directories("$ENV{NPTOOLS_DIR}/include")
link_directories("$ENV{NPTOOLS_DIR}/lib")
include_directories("$ENV{EIGEN_DIR}")


add_executable(unittest unittest.cpp linelist.cpp)
target_link_libraries(unittest gtest gtest_main)



