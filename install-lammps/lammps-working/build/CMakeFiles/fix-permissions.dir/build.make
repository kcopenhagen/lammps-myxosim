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

# Utility rule file for fix-permissions.

# Include any custom commands dependencies for this target.
include CMakeFiles/fix-permissions.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/fix-permissions.dir/progress.make

CMakeFiles/fix-permissions:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/kc32/software/lammps-working/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Fix permission errors"
	cd /home/kc32/software/lammps-working && /usr/bin/python3.6 /home/kc32/software/lammps-working/tools/coding_standard/permissions.py -f .

fix-permissions: CMakeFiles/fix-permissions
fix-permissions: CMakeFiles/fix-permissions.dir/build.make
.PHONY : fix-permissions

# Rule to build all files generated by this target.
CMakeFiles/fix-permissions.dir/build: fix-permissions
.PHONY : CMakeFiles/fix-permissions.dir/build

CMakeFiles/fix-permissions.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fix-permissions.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fix-permissions.dir/clean

CMakeFiles/fix-permissions.dir/depend:
	cd /home/kc32/software/lammps-working/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kc32/software/lammps-working/cmake /home/kc32/software/lammps-working/cmake /home/kc32/software/lammps-working/build /home/kc32/software/lammps-working/build /home/kc32/software/lammps-working/build/CMakeFiles/fix-permissions.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fix-permissions.dir/depend

