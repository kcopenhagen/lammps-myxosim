# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/kc32/software/lammps-working/cmake

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/kc32/software/lammps-working/build

# Utility rule file for pointers.h.

# Include any custom commands dependencies for this target.
include CMakeFiles/pointers.h.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/pointers.h.dir/progress.make

CMakeFiles/pointers.h: includes/lammps/pointers.h

includes/lammps/pointers.h: /home/kc32/software/lammps-working/src/pointers.h
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/kc32/software/lammps-working/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating includes/lammps/pointers.h"
	/usr/bin/cmake -E copy_if_different /home/kc32/software/lammps-working/src/pointers.h /home/kc32/software/lammps-working/build/includes/lammps/pointers.h

pointers.h: CMakeFiles/pointers.h
pointers.h: includes/lammps/pointers.h
pointers.h: CMakeFiles/pointers.h.dir/build.make
.PHONY : pointers.h

# Rule to build all files generated by this target.
CMakeFiles/pointers.h.dir/build: pointers.h
.PHONY : CMakeFiles/pointers.h.dir/build

CMakeFiles/pointers.h.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/pointers.h.dir/cmake_clean.cmake
.PHONY : CMakeFiles/pointers.h.dir/clean

CMakeFiles/pointers.h.dir/depend:
	cd /home/kc32/software/lammps-working/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kc32/software/lammps-working/cmake /home/kc32/software/lammps-working/cmake /home/kc32/software/lammps-working/build /home/kc32/software/lammps-working/build /home/kc32/software/lammps-working/build/CMakeFiles/pointers.h.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/pointers.h.dir/depend

