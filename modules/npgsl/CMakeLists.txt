# We don't do anything fancy here
cmake_minimum_required(VERSION 2.6)
set(CMAKE_BUILD_TYPE Release)

# Give your project a name
PROJECT (npgsl)

# Includes 
include_directories("$ENV{NPTOOLS_DIR}/include")
set(CMAKE_INSTALL_PREFIX "$ENV{NPTOOLS_DIR}")

# Includes
include_directories("$ENV{BOOST_DIR}/include")

# Includes
include_directories("$ENV{GSL_DIR}/include")
set(GSL_LIBS gsl.a gslcblas.a)
link_directories("$ENV{GSL_DIR}/lib")


# EXAMPLES --- YOU NEED TO CHANGE THIS AS NECESSARY
set(SOURCES gslhist.cpp)

# LIBRARIES -- an example of how to build both static and shared libs
add_library(npgsl ${SOURCES})
target_link_libraries(npgsl ${GSL_LIBS})

INSTALL (TARGETS npgsl
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
