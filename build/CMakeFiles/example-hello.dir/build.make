# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/user/Projects/Twisted_Edwards/libakrypt-0.x

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/user/Projects/Twisted_Edwards/build

# Include any dependencies generated for this target.
include CMakeFiles/example-hello.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/example-hello.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/example-hello.dir/flags.make

CMakeFiles/example-hello.dir/examples/example-hello.c.o: CMakeFiles/example-hello.dir/flags.make
CMakeFiles/example-hello.dir/examples/example-hello.c.o: /home/user/Projects/Twisted_Edwards/libakrypt-0.x/examples/example-hello.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/Projects/Twisted_Edwards/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/example-hello.dir/examples/example-hello.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/example-hello.dir/examples/example-hello.c.o   -c /home/user/Projects/Twisted_Edwards/libakrypt-0.x/examples/example-hello.c

CMakeFiles/example-hello.dir/examples/example-hello.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/example-hello.dir/examples/example-hello.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/user/Projects/Twisted_Edwards/libakrypt-0.x/examples/example-hello.c > CMakeFiles/example-hello.dir/examples/example-hello.c.i

CMakeFiles/example-hello.dir/examples/example-hello.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/example-hello.dir/examples/example-hello.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/user/Projects/Twisted_Edwards/libakrypt-0.x/examples/example-hello.c -o CMakeFiles/example-hello.dir/examples/example-hello.c.s

# Object files for target example-hello
example__hello_OBJECTS = \
"CMakeFiles/example-hello.dir/examples/example-hello.c.o"

# External object files for target example-hello
example__hello_EXTERNAL_OBJECTS =

example-hello: CMakeFiles/example-hello.dir/examples/example-hello.c.o
example-hello: CMakeFiles/example-hello.dir/build.make
example-hello: libakrypt-static.a
example-hello: CMakeFiles/example-hello.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/user/Projects/Twisted_Edwards/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable example-hello"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/example-hello.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/example-hello.dir/build: example-hello

.PHONY : CMakeFiles/example-hello.dir/build

CMakeFiles/example-hello.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/example-hello.dir/cmake_clean.cmake
.PHONY : CMakeFiles/example-hello.dir/clean

CMakeFiles/example-hello.dir/depend:
	cd /home/user/Projects/Twisted_Edwards/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/user/Projects/Twisted_Edwards/libakrypt-0.x /home/user/Projects/Twisted_Edwards/libakrypt-0.x /home/user/Projects/Twisted_Edwards/build /home/user/Projects/Twisted_Edwards/build /home/user/Projects/Twisted_Edwards/build/CMakeFiles/example-hello.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/example-hello.dir/depend

