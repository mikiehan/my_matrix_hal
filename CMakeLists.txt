cmake_minimum_required(VERSION 2.8)

# Add -Wall and -Wextra. Also,
# treat C/C++ warnings as errors if -DADM_FATAL_WARNINGS=ON.
include (cmake/FatalWarnings.cmake)
ADM_EXTRA_WARNINGS()

add_subdirectory(demos)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -I/usr/local/include/eigen-3.3.7")

