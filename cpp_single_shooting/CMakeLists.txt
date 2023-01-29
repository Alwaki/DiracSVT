cmake_minimum_required(VERSION 3.10)

project(dirac_shooting
	VERSION 1.0.0
	DESCRIPTION "Numerical solver of the dirac equation to extract binding energy between states"
	LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)
find_package(Eigen3 REQUIRED)

add_executable(${PROJECT_NAME} main.cpp)

include_directories(${Python3_INCLUDE_DIRS} ${Eigen_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${Eigen_LIBRARIES} ${PYTHON3_LIBRARIES})