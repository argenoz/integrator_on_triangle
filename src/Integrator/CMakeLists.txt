add_library(integrator integrator.cpp)

target_include_directories(integrator PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/thirdparty
)

target_link_libraries(integrator PUBLIC Fraction)