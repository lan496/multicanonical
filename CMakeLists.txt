cmake_minimum_required(VERSION 2.8)
project(multicanonical CXX)

set(CMAKE_VERBOSE_MAKEFILE TRUE)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall -Wextra")

include_directories(include)
add_subdirectory(lib)
add_executable(
  ${PROJECT_NAME}
  main.cpp
  )
target_link_libraries(${PROJECT_NAME} util)
