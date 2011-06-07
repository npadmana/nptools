# We don't do anything fancy here
cmake_minimum_required(VERSION 2.6)
set(CMAKE_BUILD_TYPE Debug)

# Give your project a name
PROJECT (baofisher)

# Includes 
include_directories("$ENV{NPTOOLS_DIR}/include")
include_directories("$ENV{EIGEN_DIR}")
link_directories("$ENV{NPTOOLS_DIR}/lib")

# Includes
include_directories("$ENV{GSL_DIR}/include")
set(GSL_LIBS gsl.a gslcblas.a)
link_directories("$ENV{GSL_DIR}/lib")


add_library(baofisher SHARED bao_forecast.c forecast.cpp fomswg.cpp)
target_link_libraries(baofisher ${GSL_LIBS})

add_executable(unittest unittest.cpp fomswg.cpp)
target_link_libraries(unittest gtest gtest_main baofisher ${GSL_LIBS})



