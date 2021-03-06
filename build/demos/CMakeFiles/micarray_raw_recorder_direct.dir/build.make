# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/pi/my_matrix_hal

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/pi/my_matrix_hal/build

# Include any dependencies generated for this target.
include demos/CMakeFiles/micarray_raw_recorder_direct.dir/depend.make

# Include the progress variables for this target.
include demos/CMakeFiles/micarray_raw_recorder_direct.dir/progress.make

# Include the compile flags for this target's objects.
include demos/CMakeFiles/micarray_raw_recorder_direct.dir/flags.make

demos/CMakeFiles/micarray_raw_recorder_direct.dir/micarray_raw_recorder_direct.cpp.o: demos/CMakeFiles/micarray_raw_recorder_direct.dir/flags.make
demos/CMakeFiles/micarray_raw_recorder_direct.dir/micarray_raw_recorder_direct.cpp.o: ../demos/micarray_raw_recorder_direct.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pi/my_matrix_hal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object demos/CMakeFiles/micarray_raw_recorder_direct.dir/micarray_raw_recorder_direct.cpp.o"
	cd /home/pi/my_matrix_hal/build/demos && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/micarray_raw_recorder_direct.dir/micarray_raw_recorder_direct.cpp.o -c /home/pi/my_matrix_hal/demos/micarray_raw_recorder_direct.cpp

demos/CMakeFiles/micarray_raw_recorder_direct.dir/micarray_raw_recorder_direct.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/micarray_raw_recorder_direct.dir/micarray_raw_recorder_direct.cpp.i"
	cd /home/pi/my_matrix_hal/build/demos && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/pi/my_matrix_hal/demos/micarray_raw_recorder_direct.cpp > CMakeFiles/micarray_raw_recorder_direct.dir/micarray_raw_recorder_direct.cpp.i

demos/CMakeFiles/micarray_raw_recorder_direct.dir/micarray_raw_recorder_direct.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/micarray_raw_recorder_direct.dir/micarray_raw_recorder_direct.cpp.s"
	cd /home/pi/my_matrix_hal/build/demos && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/pi/my_matrix_hal/demos/micarray_raw_recorder_direct.cpp -o CMakeFiles/micarray_raw_recorder_direct.dir/micarray_raw_recorder_direct.cpp.s

# Object files for target micarray_raw_recorder_direct
micarray_raw_recorder_direct_OBJECTS = \
"CMakeFiles/micarray_raw_recorder_direct.dir/micarray_raw_recorder_direct.cpp.o"

# External object files for target micarray_raw_recorder_direct
micarray_raw_recorder_direct_EXTERNAL_OBJECTS =

demos/micarray_raw_recorder_direct: demos/CMakeFiles/micarray_raw_recorder_direct.dir/micarray_raw_recorder_direct.cpp.o
demos/micarray_raw_recorder_direct: demos/CMakeFiles/micarray_raw_recorder_direct.dir/build.make
demos/micarray_raw_recorder_direct: demos/driver/libmatrix_creator_hal.so
demos/micarray_raw_recorder_direct: /usr/lib/libwiringPi.so
demos/micarray_raw_recorder_direct: /usr/lib/libwiringPiDev.so
demos/micarray_raw_recorder_direct: /usr/lib/arm-linux-gnueabihf/libcrypt.so
demos/micarray_raw_recorder_direct: /usr/lib/arm-linux-gnueabihf/libgflags.so
demos/micarray_raw_recorder_direct: /usr/lib/arm-linux-gnueabihf/libfftw3f.so
demos/micarray_raw_recorder_direct: /usr/lib/libwiringPi.so
demos/micarray_raw_recorder_direct: /usr/lib/libwiringPiDev.so
demos/micarray_raw_recorder_direct: /usr/lib/arm-linux-gnueabihf/libcrypt.so
demos/micarray_raw_recorder_direct: demos/CMakeFiles/micarray_raw_recorder_direct.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pi/my_matrix_hal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable micarray_raw_recorder_direct"
	cd /home/pi/my_matrix_hal/build/demos && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/micarray_raw_recorder_direct.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
demos/CMakeFiles/micarray_raw_recorder_direct.dir/build: demos/micarray_raw_recorder_direct

.PHONY : demos/CMakeFiles/micarray_raw_recorder_direct.dir/build

demos/CMakeFiles/micarray_raw_recorder_direct.dir/clean:
	cd /home/pi/my_matrix_hal/build/demos && $(CMAKE_COMMAND) -P CMakeFiles/micarray_raw_recorder_direct.dir/cmake_clean.cmake
.PHONY : demos/CMakeFiles/micarray_raw_recorder_direct.dir/clean

demos/CMakeFiles/micarray_raw_recorder_direct.dir/depend:
	cd /home/pi/my_matrix_hal/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pi/my_matrix_hal /home/pi/my_matrix_hal/demos /home/pi/my_matrix_hal/build /home/pi/my_matrix_hal/build/demos /home/pi/my_matrix_hal/build/demos/CMakeFiles/micarray_raw_recorder_direct.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : demos/CMakeFiles/micarray_raw_recorder_direct.dir/depend

