# We don't do anything fancy here
cmake_minimum_required(VERSION 2.6)
set(CMAKE_BUILD_TYPE Release)

# Give your project a name
PROJECT (velocity_divergence)

# Includes 
#include_directories("$ENV{NPTOOLS_DIR}/include")

# Includes
include_directories("$ENV{BOOST_DIR}/include")
include_directories("$ENV{EIGEN_DIR}")

# Includes
#include_directories("$ENV{GSL_DIR}/include")
#set(GSL_LIBS gsl gslcblas)  
#link_directories("$ENV{GSL_DIR}/lib")

include_directories("$ENV{FFTW_DIR}/include")
link_directories("$ENV{FFTW_DIR}/lib")
set(FFTW_LIBS drfftw dfftw m)


# EXAMPLES --- YOU NEED TO CHANGE THIS AS NECESSARY
add_executable (veldiv veldiv.cpp)
#add_executable (test_array test_array.cpp)
target_link_libraries(veldiv ${FFTW_LIBS})
