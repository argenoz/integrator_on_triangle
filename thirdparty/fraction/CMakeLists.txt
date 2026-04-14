add_library(Fraction Fraction.cpp)

target_include_directories(Fraction PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)