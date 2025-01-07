cmake_minimum_required(VERSION 3.5)

# Project name and optional version
project(emme VERSION 1.0.0)

set(CMAKE_CXX_FLAGS "-DEMME_EXPRESSION_TEMPLATE -DMULTI_THREAD")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -DEMME_DEBUG -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall")

# Find all source files recursively in the "src" directory
file(GLOB_RECURSE SOURCE_FILES CONFIGURE_DEPENDS "src/*.cpp")

set(INCLUDE_DIR ./include)

set(CMAKE_CXX_STANDARD 20)

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

# Add executable target
add_executable(emme ${SOURCE_FILES})

if(UNIX AND (CMAKE_GENERATOR STREQUAL "Unix Makefiles"))
    target_compile_definitions(emme PRIVATE [=[`${CMAKE_SOURCE_DIR}/build_info.sh hash`]=])
    target_compile_definitions(emme PRIVATE [=[`${CMAKE_SOURCE_DIR}/build_info.sh time`]=])
endif()

if(EMME_MKL)
    target_compile_definitions(emme PRIVATE EMME_MKL)
endif()

target_include_directories(emme PUBLIC ${INCLUDE_DIR})


# Link executable with BLAS and LAPACK
target_link_libraries(emme PUBLIC BLAS::BLAS LAPACK::LAPACK)
