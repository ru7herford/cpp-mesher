cmake_minimum_required(VERSION 3.20)
project(cpp_mesher LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_library(mesher_core STATIC
    src/predicates.cpp
    src/cdt2d.cpp
    src/refine.cpp
    src/quality.cpp
    src/io.cpp
)
target_include_directories(mesher_core PUBLIC src)
# The exact predicates rely on plain IEEE binary64 rounding; a fused
# multiply-add would silently break their error-free transforms.
target_compile_options(mesher_core PUBLIC -Wall -Wextra -ffp-contract=off)

add_executable(mesher src/main.cpp)
target_link_libraries(mesher PRIVATE mesher_core)

enable_testing()
foreach(name examples predicates fuzz)
    add_executable(test_${name} tests/test_${name}.cpp)
    target_link_libraries(test_${name} PRIVATE mesher_core)
endforeach()
add_test(NAME examples COMMAND test_examples ${CMAKE_CURRENT_SOURCE_DIR}/examples ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME predicates COMMAND test_predicates)
add_test(NAME fuzz COMMAND test_fuzz)
