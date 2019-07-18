cmake_minimum_required(VERSION 3.10.2)

set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_C_FLAGS_DEBUG "-g -O0")
set(CMAKE_C_FLAGS_RELEASE "-O3")

project(wprojection)
add_executable(wprojection main.c w_projection.c window.c utility.c fft.c complex.c)
target_link_libraries(wprojection m)