# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.19

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.19.6/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.19.6/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/andrew/Documents/CGRT/FIreInCube

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/andrew/Documents/CGRT/FIreInCube/build

# Include any dependencies generated for this target.
include CMakeFiles/fireincube.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/fireincube.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/fireincube.dir/flags.make

CMakeFiles/fireincube.dir/main.cpp.o: CMakeFiles/fireincube.dir/flags.make
CMakeFiles/fireincube.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/andrew/Documents/CGRT/FIreInCube/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/fireincube.dir/main.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fireincube.dir/main.cpp.o -c /Users/andrew/Documents/CGRT/FIreInCube/main.cpp

CMakeFiles/fireincube.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fireincube.dir/main.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/andrew/Documents/CGRT/FIreInCube/main.cpp > CMakeFiles/fireincube.dir/main.cpp.i

CMakeFiles/fireincube.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fireincube.dir/main.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/andrew/Documents/CGRT/FIreInCube/main.cpp -o CMakeFiles/fireincube.dir/main.cpp.s

# Object files for target fireincube
fireincube_OBJECTS = \
"CMakeFiles/fireincube.dir/main.cpp.o"

# External object files for target fireincube
fireincube_EXTERNAL_OBJECTS =

fireincube: CMakeFiles/fireincube.dir/main.cpp.o
fireincube: CMakeFiles/fireincube.dir/build.make
fireincube: CMakeFiles/fireincube.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/andrew/Documents/CGRT/FIreInCube/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable fireincube"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fireincube.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/fireincube.dir/build: fireincube

.PHONY : CMakeFiles/fireincube.dir/build

CMakeFiles/fireincube.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fireincube.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fireincube.dir/clean

CMakeFiles/fireincube.dir/depend:
	cd /Users/andrew/Documents/CGRT/FIreInCube/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/andrew/Documents/CGRT/FIreInCube /Users/andrew/Documents/CGRT/FIreInCube /Users/andrew/Documents/CGRT/FIreInCube/build /Users/andrew/Documents/CGRT/FIreInCube/build /Users/andrew/Documents/CGRT/FIreInCube/build/CMakeFiles/fireincube.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fireincube.dir/depend
