cmake_minimum_required(VERSION 3.28)
project(Galerkin)

add_subdirectory(src/Integrator)
add_subdirectory(thirdparty/fraction)

add_executable(Galerkin src/Galerkin.cpp)


target_link_libraries(Galerkin
    PUBLIC 
		integrator
		Fraction
)

target_include_directories(Galerkin 
	PRIVATE 
	${PROJECT_SOURCE_DIR}/include 
	${PROJECT_SOURCE_DIR}/thirdparty
)