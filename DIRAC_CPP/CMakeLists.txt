cmake_minimum_required(VERSION 3.2)

project(DIRAC_CPP)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -O3")

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

set(DIRAC_CPP_HEADERS_DIR ${PROJECT_SOURCE_DIR}/include)
set(DIRAC_CPP_BIN_DIR ${PROJECT_SOURCE_DIR}/bin)
set(DIRAC_CPP_LIB_DIR ${PROJECT_SOURCE_DIR}/lib)

include_directories(${DIRAC_CPP_LIB_DIR})
include_directories(${DIRAC_CPP_HEADERS_DIR})

add_subdirectory(src)
add_subdirectory(bin)