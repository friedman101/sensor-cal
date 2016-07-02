cmake_minimum_required(VERSION 2.8.9)
project(sensor-cal)
 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

file(GLOB SOURCES "src/*.c")

add_executable(sensor-cal ${SOURCES})
target_link_libraries(sensor-cal m gsl gslcblas)
SET(CMAKE_C_FLAGS "-std=gnu99 -Wall -pedantic -Wextra -Werror -O3")
