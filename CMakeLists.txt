cmake_minimum_required(VERSION 2.8.9)
project(sensorcal)
 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_library(sensorcal ${CMAKE_BINARY_DIR}/src/sensor_cal.c  ${CMAKE_BINARY_DIR}/src/sensor_cal.h)
target_link_libraries(sensorcal gsl gslcblas)
add_executable(test-cal ${CMAKE_BINARY_DIR}/src/test.c)
target_link_libraries(test-cal sensorcal)

SET(CMAKE_C_FLAGS "-std=gnu99 -Wall -pedantic -Wextra -Werror -O3")
