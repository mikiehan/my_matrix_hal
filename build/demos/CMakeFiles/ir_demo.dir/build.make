# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

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
CMAKE_SOURCE_DIR = /home/pi/matrix-creator-hal

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/pi/matrix-creator-hal/build

# Include any dependencies generated for this target.
include demos/CMakeFiles/ir_demo.dir/depend.make

# Include the progress variables for this target.
include demos/CMakeFiles/ir_demo.dir/progress.make

# Include the compile flags for this target's objects.
include demos/CMakeFiles/ir_demo.dir/flags.make

demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.o: demos/CMakeFiles/ir_demo.dir/flags.make
demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.o: ../demos/ir_demo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pi/matrix-creator-hal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.o"
	cd /home/pi/matrix-creator-hal/build/demos && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ir_demo.dir/ir_demo.cpp.o -c /home/pi/matrix-creator-hal/demos/ir_demo.cpp

demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ir_demo.dir/ir_demo.cpp.i"
	cd /home/pi/matrix-creator-hal/build/demos && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/pi/matrix-creator-hal/demos/ir_demo.cpp > CMakeFiles/ir_demo.dir/ir_demo.cpp.i

demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ir_demo.dir/ir_demo.cpp.s"
	cd /home/pi/matrix-creator-hal/build/demos && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/pi/matrix-creator-hal/demos/ir_demo.cpp -o CMakeFiles/ir_demo.dir/ir_demo.cpp.s

demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.o.requires:

.PHONY : demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.o.requires

demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.o.provides: demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.o.requires
	$(MAKE) -f demos/CMakeFiles/ir_demo.dir/build.make demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.o.provides.build
.PHONY : demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.o.provides

demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.o.provides.build: demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.o


# Object files for target ir_demo
ir_demo_OBJECTS = \
"CMakeFiles/ir_demo.dir/ir_demo.cpp.o"

# External object files for target ir_demo
ir_demo_EXTERNAL_OBJECTS =

demos/ir_demo: demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.o
demos/ir_demo: demos/CMakeFiles/ir_demo.dir/build.make
demos/ir_demo: demos/driver/libmatrix_creator_hal.so
demos/ir_demo: /usr/lib/libwiringPi.so
demos/ir_demo: /usr/lib/libwiringPiDev.so
demos/ir_demo: /usr/lib/arm-linux-gnueabihf/libcrypt.so
demos/ir_demo: /usr/lib/arm-linux-gnueabihf/libfftw3f.so
demos/ir_demo: /usr/lib/libwiringPi.so
demos/ir_demo: /usr/lib/libwiringPiDev.so
demos/ir_demo: /usr/lib/arm-linux-gnueabihf/libcrypt.so
demos/ir_demo: demos/CMakeFiles/ir_demo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pi/matrix-creator-hal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ir_demo"
	cd /home/pi/matrix-creator-hal/build/demos && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ir_demo.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
demos/CMakeFiles/ir_demo.dir/build: demos/ir_demo

.PHONY : demos/CMakeFiles/ir_demo.dir/build

demos/CMakeFiles/ir_demo.dir/requires: demos/CMakeFiles/ir_demo.dir/ir_demo.cpp.o.requires

.PHONY : demos/CMakeFiles/ir_demo.dir/requires

demos/CMakeFiles/ir_demo.dir/clean:
	cd /home/pi/matrix-creator-hal/build/demos && $(CMAKE_COMMAND) -P CMakeFiles/ir_demo.dir/cmake_clean.cmake
.PHONY : demos/CMakeFiles/ir_demo.dir/clean

demos/CMakeFiles/ir_demo.dir/depend:
	cd /home/pi/matrix-creator-hal/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pi/matrix-creator-hal /home/pi/matrix-creator-hal/demos /home/pi/matrix-creator-hal/build /home/pi/matrix-creator-hal/build/demos /home/pi/matrix-creator-hal/build/demos/CMakeFiles/ir_demo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : demos/CMakeFiles/ir_demo.dir/depend

