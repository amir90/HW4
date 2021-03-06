# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_SOURCE_DIR = "/home/amir/Documents/CS Masters/Robotics/HW4"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/amir/Documents/CS Masters/Robotics/HW4"

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/cmake-gui -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start "/home/amir/Documents/CS Masters/Robotics/HW4/CMakeFiles" "/home/amir/Documents/CS Masters/Robotics/HW4/CMakeFiles/progress.marks"
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start "/home/amir/Documents/CS Masters/Robotics/HW4/CMakeFiles" 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named RodPathFinder

# Build rule for target.
RodPathFinder: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 RodPathFinder
.PHONY : RodPathFinder

# fast build rule for target.
RodPathFinder/fast:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/build
.PHONY : RodPathFinder/fast

#=============================================================================
# Target rules for targets named QueryHandler

# Build rule for target.
QueryHandler: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 QueryHandler
.PHONY : QueryHandler

# fast build rule for target.
QueryHandler/fast:
	$(MAKE) -f CMakeFiles/QueryHandler.dir/build.make CMakeFiles/QueryHandler.dir/build
.PHONY : QueryHandler/fast

MyQueryHandler.o: MyQueryHandler.cpp.o

.PHONY : MyQueryHandler.o

# target to build an object file
MyQueryHandler.cpp.o:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/MyQueryHandler.cpp.o
	$(MAKE) -f CMakeFiles/QueryHandler.dir/build.make CMakeFiles/QueryHandler.dir/MyQueryHandler.cpp.o
.PHONY : MyQueryHandler.cpp.o

MyQueryHandler.i: MyQueryHandler.cpp.i

.PHONY : MyQueryHandler.i

# target to preprocess a source file
MyQueryHandler.cpp.i:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/MyQueryHandler.cpp.i
	$(MAKE) -f CMakeFiles/QueryHandler.dir/build.make CMakeFiles/QueryHandler.dir/MyQueryHandler.cpp.i
.PHONY : MyQueryHandler.cpp.i

MyQueryHandler.s: MyQueryHandler.cpp.s

.PHONY : MyQueryHandler.s

# target to generate assembly for a file
MyQueryHandler.cpp.s:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/MyQueryHandler.cpp.s
	$(MAKE) -f CMakeFiles/QueryHandler.dir/build.make CMakeFiles/QueryHandler.dir/MyQueryHandler.cpp.s
.PHONY : MyQueryHandler.cpp.s

MyRodPathFinder.o: MyRodPathFinder.cpp.o

.PHONY : MyRodPathFinder.o

# target to build an object file
MyRodPathFinder.cpp.o:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/MyRodPathFinder.cpp.o
.PHONY : MyRodPathFinder.cpp.o

MyRodPathFinder.i: MyRodPathFinder.cpp.i

.PHONY : MyRodPathFinder.i

# target to preprocess a source file
MyRodPathFinder.cpp.i:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/MyRodPathFinder.cpp.i
.PHONY : MyRodPathFinder.cpp.i

MyRodPathFinder.s: MyRodPathFinder.cpp.s

.PHONY : MyRodPathFinder.s

# target to generate assembly for a file
MyRodPathFinder.cpp.s:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/MyRodPathFinder.cpp.s
.PHONY : MyRodPathFinder.cpp.s

NaiveQueryHandler.o: NaiveQueryHandler.cpp.o

.PHONY : NaiveQueryHandler.o

# target to build an object file
NaiveQueryHandler.cpp.o:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/NaiveQueryHandler.cpp.o
	$(MAKE) -f CMakeFiles/QueryHandler.dir/build.make CMakeFiles/QueryHandler.dir/NaiveQueryHandler.cpp.o
.PHONY : NaiveQueryHandler.cpp.o

NaiveQueryHandler.i: NaiveQueryHandler.cpp.i

.PHONY : NaiveQueryHandler.i

# target to preprocess a source file
NaiveQueryHandler.cpp.i:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/NaiveQueryHandler.cpp.i
	$(MAKE) -f CMakeFiles/QueryHandler.dir/build.make CMakeFiles/QueryHandler.dir/NaiveQueryHandler.cpp.i
.PHONY : NaiveQueryHandler.cpp.i

NaiveQueryHandler.s: NaiveQueryHandler.cpp.s

.PHONY : NaiveQueryHandler.s

# target to generate assembly for a file
NaiveQueryHandler.cpp.s:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/NaiveQueryHandler.cpp.s
	$(MAKE) -f CMakeFiles/QueryHandler.dir/build.make CMakeFiles/QueryHandler.dir/NaiveQueryHandler.cpp.s
.PHONY : NaiveQueryHandler.cpp.s

Path.o: Path.cpp.o

.PHONY : Path.o

# target to build an object file
Path.cpp.o:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/Path.cpp.o
.PHONY : Path.cpp.o

Path.i: Path.cpp.i

.PHONY : Path.i

# target to preprocess a source file
Path.cpp.i:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/Path.cpp.i
.PHONY : Path.cpp.i

Path.s: Path.cpp.s

.PHONY : Path.s

# target to generate assembly for a file
Path.cpp.s:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/Path.cpp.s
.PHONY : Path.cpp.s

QueryHandler.o: QueryHandler.cpp.o

.PHONY : QueryHandler.o

# target to build an object file
QueryHandler.cpp.o:
	$(MAKE) -f CMakeFiles/QueryHandler.dir/build.make CMakeFiles/QueryHandler.dir/QueryHandler.cpp.o
.PHONY : QueryHandler.cpp.o

QueryHandler.i: QueryHandler.cpp.i

.PHONY : QueryHandler.i

# target to preprocess a source file
QueryHandler.cpp.i:
	$(MAKE) -f CMakeFiles/QueryHandler.dir/build.make CMakeFiles/QueryHandler.dir/QueryHandler.cpp.i
.PHONY : QueryHandler.cpp.i

QueryHandler.s: QueryHandler.cpp.s

.PHONY : QueryHandler.s

# target to generate assembly for a file
QueryHandler.cpp.s:
	$(MAKE) -f CMakeFiles/QueryHandler.dir/build.make CMakeFiles/QueryHandler.dir/QueryHandler.cpp.s
.PHONY : QueryHandler.cpp.s

RodPathFinder.o: RodPathFinder.cpp.o

.PHONY : RodPathFinder.o

# target to build an object file
RodPathFinder.cpp.o:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/RodPathFinder.cpp.o
.PHONY : RodPathFinder.cpp.o

RodPathFinder.i: RodPathFinder.cpp.i

.PHONY : RodPathFinder.i

# target to preprocess a source file
RodPathFinder.cpp.i:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/RodPathFinder.cpp.i
.PHONY : RodPathFinder.cpp.i

RodPathFinder.s: RodPathFinder.cpp.s

.PHONY : RodPathFinder.s

# target to generate assembly for a file
RodPathFinder.cpp.s:
	$(MAKE) -f CMakeFiles/RodPathFinder.dir/build.make CMakeFiles/RodPathFinder.dir/RodPathFinder.cpp.s
.PHONY : RodPathFinder.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... RodPathFinder"
	@echo "... QueryHandler"
	@echo "... MyQueryHandler.o"
	@echo "... MyQueryHandler.i"
	@echo "... MyQueryHandler.s"
	@echo "... MyRodPathFinder.o"
	@echo "... MyRodPathFinder.i"
	@echo "... MyRodPathFinder.s"
	@echo "... NaiveQueryHandler.o"
	@echo "... NaiveQueryHandler.i"
	@echo "... NaiveQueryHandler.s"
	@echo "... Path.o"
	@echo "... Path.i"
	@echo "... Path.s"
	@echo "... QueryHandler.o"
	@echo "... QueryHandler.i"
	@echo "... QueryHandler.s"
	@echo "... RodPathFinder.o"
	@echo "... RodPathFinder.i"
	@echo "... RodPathFinder.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

