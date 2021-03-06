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
include demos/CMakeFiles/mic_demo_fir_direct.dir/depend.make

# Include the progress variables for this target.
include demos/CMakeFiles/mic_demo_fir_direct.dir/progress.make

# Include the compile flags for this target's objects.
include demos/CMakeFiles/mic_demo_fir_direct.dir/flags.make

demos/CMakeFiles/mic_demo_fir_direct.dir/mic_demo_fir_direct.cpp.o: demos/CMakeFiles/mic_demo_fir_direct.dir/flags.make
demos/CMakeFiles/mic_demo_fir_direct.dir/mic_demo_fir_direct.cpp.o: ../demos/mic_demo_fir_direct.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pi/my_matrix_hal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object demos/CMakeFiles/mic_demo_fir_direct.dir/mic_demo_fir_direct.cpp.o"
	cd /home/pi/my_matrix_hal/build/demos && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mic_demo_fir_direct.dir/mic_demo_fir_direct.cpp.o -c /home/pi/my_matrix_hal/demos/mic_demo_fir_direct.cpp

demos/CMakeFiles/mic_demo_fir_direct.dir/mic_demo_fir_direct.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mic_demo_fir_direct.dir/mic_demo_fir_direct.cpp.i"
	cd /home/pi/my_matrix_hal/build/demos && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/pi/my_matrix_hal/demos/mic_demo_fir_direct.cpp > CMakeFiles/mic_demo_fir_direct.dir/mic_demo_fir_direct.cpp.i

demos/CMakeFiles/mic_demo_fir_direct.dir/mic_demo_fir_direct.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mic_demo_fir_direct.dir/mic_demo_fir_direct.cpp.s"
	cd /home/pi/my_matrix_hal/build/demos && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/pi/my_matrix_hal/demos/mic_demo_fir_direct.cpp -o CMakeFiles/mic_demo_fir_direct.dir/mic_demo_fir_direct.cpp.s

# Object files for target mic_demo_fir_direct
mic_demo_fir_direct_OBJECTS = \
"CMakeFiles/mic_demo_fir_direct.dir/mic_demo_fir_direct.cpp.o"

# External object files for target mic_demo_fir_direct
mic_demo_fir_direct_EXTERNAL_OBJECTS =

demos/mic_demo_fir_direct: demos/CMakeFiles/mic_demo_fir_direct.dir/mic_demo_fir_direct.cpp.o
demos/mic_demo_fir_direct: demos/CMakeFiles/mic_demo_fir_direct.dir/build.make
demos/mic_demo_fir_direct: demos/driver/libmatrix_creator_hal.so
demos/mic_demo_fir_direct: /usr/lib/libwiringPi.so
demos/mic_demo_fir_direct: /usr/lib/libwiringPiDev.so
demos/mic_demo_fir_direct: /usr/lib/arm-linux-gnueabihf/libcrypt.so
demos/mic_demo_fir_direct: /usr/lib/arm-linux-gnueabihf/libgflags.so
demos/mic_demo_fir_direct: /usr/lib/arm-linux-gnueabihf/libfftw3f.so
demos/mic_demo_fir_direct: /usr/lib/libwiringPi.so
demos/mic_demo_fir_direct: /usr/lib/libwiringPiDev.so
demos/mic_demo_fir_direct: /usr/lib/arm-linux-gnueabihf/libcrypt.so
demos/mic_demo_fir_direct: demos/CMakeFiles/mic_demo_fir_direct.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pi/my_matrix_hal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable mic_demo_fir_direct"
	cd /home/pi/my_matrix_hal/build/demos && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mic_demo_fir_direct.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
demos/CMakeFiles/mic_demo_fir_direct.dir/build: demos/mic_demo_fir_direct

.PHONY : demos/CMakeFiles/mic_demo_fir_direct.dir/build

demos/CMakeFiles/mic_demo_fir_direct.dir/clean:
	cd /home/pi/my_matrix_hal/build/demos && $(CMAKE_COMMAND) -P CMakeFiles/mic_demo_fir_direct.dir/cmake_clean.cmake
.PHONY : demos/CMakeFiles/mic_demo_fir_direct.dir/clean

demos/CMakeFiles/mic_demo_fir_direct.dir/depend:
	cd /home/pi/my_matrix_hal/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pi/my_matrix_hal /home/pi/my_matrix_hal/demos /home/pi/my_matrix_hal/build /home/pi/my_matrix_hal/build/demos /home/pi/my_matrix_hal/build/demos/CMakeFiles/mic_demo_fir_direct.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : demos/CMakeFiles/mic_demo_fir_direct.dir/depend

