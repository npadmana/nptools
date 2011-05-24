# We don't do anything fancy here
cmake_minimum_required(VERSION 2.6)
set(CMAKE_BUILD_TYPE Debug)

# Give your project a name
PROJECT (wfirst_fisher)

# Includes 
include_directories("$ENV{NPTOOLS_DIR}/include")

# Includes
include_directories("$ENV{BOOST_DIR}/include")

# Includes
include_directories("$ENV{GSL_DIR}/include")
set(GSL_LIBS gsl.a gslcblas.a)
link_directories("$ENV{GSL_DIR}/lib")


# EXAMPLES --- YOU NEED TO CHANGE THIS AS NECESSARY
set(SOURCES fomswg.cpp fisher_utils.cpp wfirst_detf.cpp)


# LIBRARIES -- an example of how to build both static and shared libs
add_library(wfirst_fisher ${SOURCES})
target_link_libraries(wfirst_fisher ${GSL_LIBS})

# EXECUTABLES
add_executable(fomswg_stage3 fomswg_stage3.cpp)
target_link_libraries(fomswg_stage3 wfirst_fisher)

add_executable(fomswg_to_detf fomswg_to_detf.cpp)
target_link_libraries(fomswg_to_detf wfirst_fisher)

add_executable(combine_detf_fisher combine_detf_fisher.cpp)
target_link_libraries(combine_detf_fisher wfirst_fisher)

add_executable(mk_sn_fisher mk_sn_fisher.cpp)
target_link_libraries(mk_sn_fisher wfirst_fisher)

add_executable(test_detf_cosmology test_detf_cosmology.cpp)
target_link_libraries(test_detf_cosmology wfirst_fisher)
