cmake_minimum_required(VERSION 3.0)

# Project name and optional version
project(emme VERSION 1.0.0)

set(MULTI_THREAD ON)

set(CMAKE_CXX_FLAGS "-DEMME_EXPRESSION_TEMPLATE -DMULTI_THREAD")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -DEMME_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall")

# Find all source files recursively in the "src" directory
file(GLOB_RECURSE SOURCE_FILES CONFIGURE_DEPENDS "src/*.cpp")

set(INCLUDE_DIR ./include)

set(CMAKE_CXX_STANDARD 20)

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

option(BUILD_DEBUG "Build with debug symbols")
option(BUILD_RELEASE "Build for production (optimized)")

set(MULTI_THREAD VALUE)  # Sets the TESTTTT macro to the specified VALUE

# if(BUILD_DEBUG)
#     # Add compiler flags for debug symbols
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
# endif()

# Add executable target
add_executable(emme ${SOURCE_FILES})

target_include_directories(emme PUBLIC ${INCLUDE_DIR})


# Link executable with BLAS and LAPACK
target_link_libraries(emme PUBLIC BLAS::BLAS LAPACK::LAPACK)
