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

# Utility rule file for gitversion.

# Include any custom commands dependencies for this target.
include CMakeFiles/gitversion.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/gitversion.dir/progress.make

CMakeFiles/gitversion:
	/usr/bin/cmake -DLAMMPS_DIR="/home/kc32/software/lammps-working" -DGIT_EXECUTABLE="/usr/bin/git" -DGIT_FOUND="TRUE" -DLAMMPS_STYLE_HEADERS_DIR="/home/kc32/software/lammps-working/build/styles" -P /home/kc32/software/lammps-working/cmake/Modules/generate_lmpgitversion.cmake

gitversion: CMakeFiles/gitversion
gitversion: CMakeFiles/gitversion.dir/build.make
.PHONY : gitversion

# Rule to build all files generated by this target.
CMakeFiles/gitversion.dir/build: gitversion
.PHONY : CMakeFiles/gitversion.dir/build

CMakeFiles/gitversion.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/gitversion.dir/cmake_clean.cmake
.PHONY : CMakeFiles/gitversion.dir/clean

CMakeFiles/gitversion.dir/depend:
	cd /home/kc32/software/lammps-working/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kc32/software/lammps-working/cmake /home/kc32/software/lammps-working/cmake /home/kc32/software/lammps-working/build /home/kc32/software/lammps-working/build /home/kc32/software/lammps-working/build/CMakeFiles/gitversion.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/gitversion.dir/depend

