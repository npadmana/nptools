# We don't do anything fancy here
cmake_minimum_required(VERSION 2.6)
set(CMAKE_BUILD_TYPE Debug)

# Give your project a name
PROJECT (velocity_divergence)

# Includes 
#include_directories("$ENV{NPTOOLS_DIR}/include")

# Includes
include_directories("$ENV{BOOST_DIR}/include")

# Includes
#include_directories("$ENV{GSL_DIR}/include")
#set(GSL_LIBS gsl gslcblas)  
#link_directories("$ENV{GSL_DIR}/lib")


# EXAMPLES --- YOU NEED TO CHANGE THIS AS NECESSARY
add_executable (veldiv veldiv.cpp)
add_executable (test_array test_array.cpp)
