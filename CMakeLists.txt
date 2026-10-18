cmake_minimum_required(VERSION 3.16)
project(quadrature_workbench LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_library(qw STATIC src/quadrature.cpp)
target_include_directories(qw PUBLIC src)
if(MSVC)
    target_compile_options(qw PRIVATE /W4)
else()
    target_compile_options(qw PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(test_quadrature tests/test_quadrature.cpp)
target_link_libraries(test_quadrature PRIVATE qw)
if(MSVC)
    target_compile_options(test_quadrature PRIVATE /W4)
else()
    target_compile_options(test_quadrature PRIVATE -Wall -Wextra -Wpedantic)
endif()

enable_testing()
add_test(NAME quadrature_checks COMMAND test_quadrature)
