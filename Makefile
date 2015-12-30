# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

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
CMAKE_SOURCE_DIR = /home/antoniofisk/Desktop/Uni/MasterThesis/Project

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/antoniofisk/Desktop/Uni/MasterThesis/Project

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running interactive CMake command-line interface..."
	/usr/bin/cmake -i .
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/antoniofisk/Desktop/Uni/MasterThesis/Project/CMakeFiles /home/antoniofisk/Desktop/Uni/MasterThesis/Project/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/antoniofisk/Desktop/Uni/MasterThesis/Project/CMakeFiles 0
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
# Target rules for targets named EMILI

# Build rule for target.
EMILI: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 EMILI
.PHONY : EMILI

# fast build rule for target.
EMILI/fast:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/build
.PHONY : EMILI/fast

ROADF/driver.o: ROADF/driver.cpp.o
.PHONY : ROADF/driver.o

# target to build an object file
ROADF/driver.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/driver.cpp.o
.PHONY : ROADF/driver.cpp.o

ROADF/driver.i: ROADF/driver.cpp.i
.PHONY : ROADF/driver.i

# target to preprocess a source file
ROADF/driver.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/driver.cpp.i
.PHONY : ROADF/driver.cpp.i

ROADF/driver.s: ROADF/driver.cpp.s
.PHONY : ROADF/driver.s

# target to generate assembly for a file
ROADF/driver.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/driver.cpp.s
.PHONY : ROADF/driver.cpp.s

ROADF/instance.o: ROADF/instance.cpp.o
.PHONY : ROADF/instance.o

# target to build an object file
ROADF/instance.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/instance.cpp.o
.PHONY : ROADF/instance.cpp.o

ROADF/instance.i: ROADF/instance.cpp.i
.PHONY : ROADF/instance.i

# target to preprocess a source file
ROADF/instance.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/instance.cpp.i
.PHONY : ROADF/instance.cpp.i

ROADF/instance.s: ROADF/instance.cpp.s
.PHONY : ROADF/instance.s

# target to generate assembly for a file
ROADF/instance.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/instance.cpp.s
.PHONY : ROADF/instance.cpp.s

ROADF/location.o: ROADF/location.cpp.o
.PHONY : ROADF/location.o

# target to build an object file
ROADF/location.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/location.cpp.o
.PHONY : ROADF/location.cpp.o

ROADF/location.i: ROADF/location.cpp.i
.PHONY : ROADF/location.i

# target to preprocess a source file
ROADF/location.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/location.cpp.i
.PHONY : ROADF/location.cpp.i

ROADF/location.s: ROADF/location.cpp.s
.PHONY : ROADF/location.s

# target to generate assembly for a file
ROADF/location.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/location.cpp.s
.PHONY : ROADF/location.cpp.s

ROADF/main.o: ROADF/main.cpp.o
.PHONY : ROADF/main.o

# target to build an object file
ROADF/main.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/main.cpp.o
.PHONY : ROADF/main.cpp.o

ROADF/main.i: ROADF/main.cpp.i
.PHONY : ROADF/main.i

# target to preprocess a source file
ROADF/main.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/main.cpp.i
.PHONY : ROADF/main.cpp.i

ROADF/main.s: ROADF/main.cpp.s
.PHONY : ROADF/main.s

# target to generate assembly for a file
ROADF/main.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/main.cpp.s
.PHONY : ROADF/main.cpp.s

ROADF/matrix.o: ROADF/matrix.cpp.o
.PHONY : ROADF/matrix.o

# target to build an object file
ROADF/matrix.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/matrix.cpp.o
.PHONY : ROADF/matrix.cpp.o

ROADF/matrix.i: ROADF/matrix.cpp.i
.PHONY : ROADF/matrix.i

# target to preprocess a source file
ROADF/matrix.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/matrix.cpp.i
.PHONY : ROADF/matrix.cpp.i

ROADF/matrix.s: ROADF/matrix.cpp.s
.PHONY : ROADF/matrix.s

# target to generate assembly for a file
ROADF/matrix.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/matrix.cpp.s
.PHONY : ROADF/matrix.cpp.s

ROADF/operation.o: ROADF/operation.cpp.o
.PHONY : ROADF/operation.o

# target to build an object file
ROADF/operation.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/operation.cpp.o
.PHONY : ROADF/operation.cpp.o

ROADF/operation.i: ROADF/operation.cpp.i
.PHONY : ROADF/operation.i

# target to preprocess a source file
ROADF/operation.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/operation.cpp.i
.PHONY : ROADF/operation.cpp.i

ROADF/operation.s: ROADF/operation.cpp.s
.PHONY : ROADF/operation.s

# target to generate assembly for a file
ROADF/operation.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/operation.cpp.s
.PHONY : ROADF/operation.cpp.s

ROADF/solution.o: ROADF/solution.cpp.o
.PHONY : ROADF/solution.o

# target to build an object file
ROADF/solution.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/solution.cpp.o
.PHONY : ROADF/solution.cpp.o

ROADF/solution.i: ROADF/solution.cpp.i
.PHONY : ROADF/solution.i

# target to preprocess a source file
ROADF/solution.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/solution.cpp.i
.PHONY : ROADF/solution.cpp.i

ROADF/solution.s: ROADF/solution.cpp.s
.PHONY : ROADF/solution.s

# target to generate assembly for a file
ROADF/solution.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/solution.cpp.s
.PHONY : ROADF/solution.cpp.s

ROADF/tinystr.o: ROADF/tinystr.cpp.o
.PHONY : ROADF/tinystr.o

# target to build an object file
ROADF/tinystr.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/tinystr.cpp.o
.PHONY : ROADF/tinystr.cpp.o

ROADF/tinystr.i: ROADF/tinystr.cpp.i
.PHONY : ROADF/tinystr.i

# target to preprocess a source file
ROADF/tinystr.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/tinystr.cpp.i
.PHONY : ROADF/tinystr.cpp.i

ROADF/tinystr.s: ROADF/tinystr.cpp.s
.PHONY : ROADF/tinystr.s

# target to generate assembly for a file
ROADF/tinystr.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/tinystr.cpp.s
.PHONY : ROADF/tinystr.cpp.s

ROADF/tinyxml.o: ROADF/tinyxml.cpp.o
.PHONY : ROADF/tinyxml.o

# target to build an object file
ROADF/tinyxml.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/tinyxml.cpp.o
.PHONY : ROADF/tinyxml.cpp.o

ROADF/tinyxml.i: ROADF/tinyxml.cpp.i
.PHONY : ROADF/tinyxml.i

# target to preprocess a source file
ROADF/tinyxml.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/tinyxml.cpp.i
.PHONY : ROADF/tinyxml.cpp.i

ROADF/tinyxml.s: ROADF/tinyxml.cpp.s
.PHONY : ROADF/tinyxml.s

# target to generate assembly for a file
ROADF/tinyxml.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/tinyxml.cpp.s
.PHONY : ROADF/tinyxml.cpp.s

ROADF/tinyxmlerror.o: ROADF/tinyxmlerror.cpp.o
.PHONY : ROADF/tinyxmlerror.o

# target to build an object file
ROADF/tinyxmlerror.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/tinyxmlerror.cpp.o
.PHONY : ROADF/tinyxmlerror.cpp.o

ROADF/tinyxmlerror.i: ROADF/tinyxmlerror.cpp.i
.PHONY : ROADF/tinyxmlerror.i

# target to preprocess a source file
ROADF/tinyxmlerror.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/tinyxmlerror.cpp.i
.PHONY : ROADF/tinyxmlerror.cpp.i

ROADF/tinyxmlerror.s: ROADF/tinyxmlerror.cpp.s
.PHONY : ROADF/tinyxmlerror.s

# target to generate assembly for a file
ROADF/tinyxmlerror.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/tinyxmlerror.cpp.s
.PHONY : ROADF/tinyxmlerror.cpp.s

ROADF/tinyxmlparser.o: ROADF/tinyxmlparser.cpp.o
.PHONY : ROADF/tinyxmlparser.o

# target to build an object file
ROADF/tinyxmlparser.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/tinyxmlparser.cpp.o
.PHONY : ROADF/tinyxmlparser.cpp.o

ROADF/tinyxmlparser.i: ROADF/tinyxmlparser.cpp.i
.PHONY : ROADF/tinyxmlparser.i

# target to preprocess a source file
ROADF/tinyxmlparser.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/tinyxmlparser.cpp.i
.PHONY : ROADF/tinyxmlparser.cpp.i

ROADF/tinyxmlparser.s: ROADF/tinyxmlparser.cpp.s
.PHONY : ROADF/tinyxmlparser.s

# target to generate assembly for a file
ROADF/tinyxmlparser.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/tinyxmlparser.cpp.s
.PHONY : ROADF/tinyxmlparser.cpp.s

ROADF/trailer.o: ROADF/trailer.cpp.o
.PHONY : ROADF/trailer.o

# target to build an object file
ROADF/trailer.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/trailer.cpp.o
.PHONY : ROADF/trailer.cpp.o

ROADF/trailer.i: ROADF/trailer.cpp.i
.PHONY : ROADF/trailer.i

# target to preprocess a source file
ROADF/trailer.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/trailer.cpp.i
.PHONY : ROADF/trailer.cpp.i

ROADF/trailer.s: ROADF/trailer.cpp.s
.PHONY : ROADF/trailer.s

# target to generate assembly for a file
ROADF/trailer.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/ROADF/trailer.cpp.s
.PHONY : ROADF/trailer.cpp.s

emilibase.o: emilibase.cpp.o
.PHONY : emilibase.o

# target to build an object file
emilibase.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/emilibase.cpp.o
.PHONY : emilibase.cpp.o

emilibase.i: emilibase.cpp.i
.PHONY : emilibase.i

# target to preprocess a source file
emilibase.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/emilibase.cpp.i
.PHONY : emilibase.cpp.i

emilibase.s: emilibase.cpp.s
.PHONY : emilibase.s

# target to generate assembly for a file
emilibase.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/emilibase.cpp.s
.PHONY : emilibase.cpp.s

generalParser.o: generalParser.cpp.o
.PHONY : generalParser.o

# target to build an object file
generalParser.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/generalParser.cpp.o
.PHONY : generalParser.cpp.o

generalParser.i: generalParser.cpp.i
.PHONY : generalParser.i

# target to preprocess a source file
generalParser.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/generalParser.cpp.i
.PHONY : generalParser.cpp.i

generalParser.s: generalParser.cpp.s
.PHONY : generalParser.s

# target to generate assembly for a file
generalParser.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/generalParser.cpp.s
.PHONY : generalParser.cpp.s

irp.o: irp.cpp.o
.PHONY : irp.o

# target to build an object file
irp.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/irp.cpp.o
.PHONY : irp.cpp.o

irp.i: irp.cpp.i
.PHONY : irp.i

# target to preprocess a source file
irp.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/irp.cpp.i
.PHONY : irp.cpp.i

irp.s: irp.cpp.s
.PHONY : irp.s

# target to generate assembly for a file
irp.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/irp.cpp.s
.PHONY : irp.cpp.s

irpparser.o: irpparser.cpp.o
.PHONY : irpparser.o

# target to build an object file
irpparser.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/irpparser.cpp.o
.PHONY : irpparser.cpp.o

irpparser.i: irpparser.cpp.i
.PHONY : irpparser.i

# target to preprocess a source file
irpparser.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/irpparser.cpp.i
.PHONY : irpparser.cpp.i

irpparser.s: irpparser.cpp.s
.PHONY : irpparser.s

# target to generate assembly for a file
irpparser.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/irpparser.cpp.s
.PHONY : irpparser.cpp.s

main.o: main.cpp.o
.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i
.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s
.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) -f CMakeFiles/EMILI.dir/build.make CMakeFiles/EMILI.dir/main.cpp.s
.PHONY : main.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... EMILI"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... ROADF/driver.o"
	@echo "... ROADF/driver.i"
	@echo "... ROADF/driver.s"
	@echo "... ROADF/instance.o"
	@echo "... ROADF/instance.i"
	@echo "... ROADF/instance.s"
	@echo "... ROADF/location.o"
	@echo "... ROADF/location.i"
	@echo "... ROADF/location.s"
	@echo "... ROADF/main.o"
	@echo "... ROADF/main.i"
	@echo "... ROADF/main.s"
	@echo "... ROADF/matrix.o"
	@echo "... ROADF/matrix.i"
	@echo "... ROADF/matrix.s"
	@echo "... ROADF/operation.o"
	@echo "... ROADF/operation.i"
	@echo "... ROADF/operation.s"
	@echo "... ROADF/solution.o"
	@echo "... ROADF/solution.i"
	@echo "... ROADF/solution.s"
	@echo "... ROADF/tinystr.o"
	@echo "... ROADF/tinystr.i"
	@echo "... ROADF/tinystr.s"
	@echo "... ROADF/tinyxml.o"
	@echo "... ROADF/tinyxml.i"
	@echo "... ROADF/tinyxml.s"
	@echo "... ROADF/tinyxmlerror.o"
	@echo "... ROADF/tinyxmlerror.i"
	@echo "... ROADF/tinyxmlerror.s"
	@echo "... ROADF/tinyxmlparser.o"
	@echo "... ROADF/tinyxmlparser.i"
	@echo "... ROADF/tinyxmlparser.s"
	@echo "... ROADF/trailer.o"
	@echo "... ROADF/trailer.i"
	@echo "... ROADF/trailer.s"
	@echo "... emilibase.o"
	@echo "... emilibase.i"
	@echo "... emilibase.s"
	@echo "... generalParser.o"
	@echo "... generalParser.i"
	@echo "... generalParser.s"
	@echo "... irp.o"
	@echo "... irp.i"
	@echo "... irp.s"
	@echo "... irpparser.o"
	@echo "... irpparser.i"
	@echo "... irpparser.s"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

