cmake_minimum_required(VERSION 2.8.9)
project(test-cal)
 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

file(GLOB SOURCES "src/*.c")

add_executable(test-cal ${SOURCES})
target_link_libraries(test-cal m gsl gslcblas)
SET(CMAKE_C_FLAGS "-std=gnu99 -Wall -pedantic -Wextra -Werror -O3")
