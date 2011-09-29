# We don't do anything fancy here
cmake_minimum_required(VERSION 2.6)
set(CMAKE_BUILD_TYPE Debug)

# Give your project a name
PROJECT (phoenix)

# Includes 
include_directories("$ENV{NPTOOLS_DIR}/include")
include_directories("$ENV{EIGEN_DIR}")
include_directories("$ENV{BOOST_DIR}/include")
link_directories("$ENV{NPTOOLS_DIR}/lib")

# Includes
include_directories("$ENV{GSL_DIR}/include")
set(GSL_LIBS gsl.a gslcblas.a)
link_directories("$ENV{GSL_DIR}/lib")

add_executable(test1 test1.cpp)



