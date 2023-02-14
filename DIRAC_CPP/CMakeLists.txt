cmake_minimum_required(VERSION 3.25.2)

project(DIRAC_SOLVER
	VERSION 1.0.0
	LANGUAGES CXX
)

set(default_build_type "Release")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -O3")

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

set(DIRAC_CPP_HEADERS_DIR ${PROJECT_SOURCE_DIR}/include)
set(DIRAC_CPP_BIN_DIR ${PROJECT_SOURCE_DIR}/bin)
set(DIRAC_CPP_LIB_DIR ${PROJECT_SOURCE_DIR}/lib)
set(DIRAC_CPP_SRC_DIR ${PROJECT_SOURCE_DIR}/src)

include_directories(${DIRAC_CPP_LIB_DIR})
include_directories(${DIRAC_CPP_HEADERS_DIR})
include_directories(${DIRAC_CPP_BIN_DIR})
include_directories(${DIRAC_CPP_SRC_DIR})


add_executable(${PROJECT_NAME} src/main.cpp src/solver.cpp src/util.cpp)