# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/baum/revSUPRA

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/baum/revSUPRA/build

# Include any dependencies generated for this target.
include CMakeFiles/supra.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/supra.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/supra.dir/flags.make

CMakeFiles/supra.dir/main_files/main.cpp.o: CMakeFiles/supra.dir/flags.make
CMakeFiles/supra.dir/main_files/main.cpp.o: ../main_files/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/baum/revSUPRA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/supra.dir/main_files/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/supra.dir/main_files/main.cpp.o -c /home/baum/revSUPRA/main_files/main.cpp

CMakeFiles/supra.dir/main_files/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/supra.dir/main_files/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/baum/revSUPRA/main_files/main.cpp > CMakeFiles/supra.dir/main_files/main.cpp.i

CMakeFiles/supra.dir/main_files/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/supra.dir/main_files/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/baum/revSUPRA/main_files/main.cpp -o CMakeFiles/supra.dir/main_files/main.cpp.s

CMakeFiles/supra.dir/src/analyzer.cpp.o: CMakeFiles/supra.dir/flags.make
CMakeFiles/supra.dir/src/analyzer.cpp.o: ../src/analyzer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/baum/revSUPRA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/supra.dir/src/analyzer.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/supra.dir/src/analyzer.cpp.o -c /home/baum/revSUPRA/src/analyzer.cpp

CMakeFiles/supra.dir/src/analyzer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/supra.dir/src/analyzer.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/baum/revSUPRA/src/analyzer.cpp > CMakeFiles/supra.dir/src/analyzer.cpp.i

CMakeFiles/supra.dir/src/analyzer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/supra.dir/src/analyzer.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/baum/revSUPRA/src/analyzer.cpp -o CMakeFiles/supra.dir/src/analyzer.cpp.s

CMakeFiles/supra.dir/src/conformergenerator.cpp.o: CMakeFiles/supra.dir/flags.make
CMakeFiles/supra.dir/src/conformergenerator.cpp.o: ../src/conformergenerator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/baum/revSUPRA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/supra.dir/src/conformergenerator.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/supra.dir/src/conformergenerator.cpp.o -c /home/baum/revSUPRA/src/conformergenerator.cpp

CMakeFiles/supra.dir/src/conformergenerator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/supra.dir/src/conformergenerator.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/baum/revSUPRA/src/conformergenerator.cpp > CMakeFiles/supra.dir/src/conformergenerator.cpp.i

CMakeFiles/supra.dir/src/conformergenerator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/supra.dir/src/conformergenerator.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/baum/revSUPRA/src/conformergenerator.cpp -o CMakeFiles/supra.dir/src/conformergenerator.cpp.s

CMakeFiles/supra.dir/src/helper.cpp.o: CMakeFiles/supra.dir/flags.make
CMakeFiles/supra.dir/src/helper.cpp.o: ../src/helper.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/baum/revSUPRA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/supra.dir/src/helper.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/supra.dir/src/helper.cpp.o -c /home/baum/revSUPRA/src/helper.cpp

CMakeFiles/supra.dir/src/helper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/supra.dir/src/helper.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/baum/revSUPRA/src/helper.cpp > CMakeFiles/supra.dir/src/helper.cpp.i

CMakeFiles/supra.dir/src/helper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/supra.dir/src/helper.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/baum/revSUPRA/src/helper.cpp -o CMakeFiles/supra.dir/src/helper.cpp.s

CMakeFiles/supra.dir/src/hungarian.cpp.o: CMakeFiles/supra.dir/flags.make
CMakeFiles/supra.dir/src/hungarian.cpp.o: ../src/hungarian.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/baum/revSUPRA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/supra.dir/src/hungarian.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/supra.dir/src/hungarian.cpp.o -c /home/baum/revSUPRA/src/hungarian.cpp

CMakeFiles/supra.dir/src/hungarian.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/supra.dir/src/hungarian.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/baum/revSUPRA/src/hungarian.cpp > CMakeFiles/supra.dir/src/hungarian.cpp.i

CMakeFiles/supra.dir/src/hungarian.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/supra.dir/src/hungarian.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/baum/revSUPRA/src/hungarian.cpp -o CMakeFiles/supra.dir/src/hungarian.cpp.s

CMakeFiles/supra.dir/src/structure.cpp.o: CMakeFiles/supra.dir/flags.make
CMakeFiles/supra.dir/src/structure.cpp.o: ../src/structure.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/baum/revSUPRA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/supra.dir/src/structure.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/supra.dir/src/structure.cpp.o -c /home/baum/revSUPRA/src/structure.cpp

CMakeFiles/supra.dir/src/structure.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/supra.dir/src/structure.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/baum/revSUPRA/src/structure.cpp > CMakeFiles/supra.dir/src/structure.cpp.i

CMakeFiles/supra.dir/src/structure.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/supra.dir/src/structure.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/baum/revSUPRA/src/structure.cpp -o CMakeFiles/supra.dir/src/structure.cpp.s

CMakeFiles/supra.dir/src/symmetry.cpp.o: CMakeFiles/supra.dir/flags.make
CMakeFiles/supra.dir/src/symmetry.cpp.o: ../src/symmetry.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/baum/revSUPRA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/supra.dir/src/symmetry.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/supra.dir/src/symmetry.cpp.o -c /home/baum/revSUPRA/src/symmetry.cpp

CMakeFiles/supra.dir/src/symmetry.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/supra.dir/src/symmetry.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/baum/revSUPRA/src/symmetry.cpp > CMakeFiles/supra.dir/src/symmetry.cpp.i

CMakeFiles/supra.dir/src/symmetry.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/supra.dir/src/symmetry.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/baum/revSUPRA/src/symmetry.cpp -o CMakeFiles/supra.dir/src/symmetry.cpp.s

# Object files for target supra
supra_OBJECTS = \
"CMakeFiles/supra.dir/main_files/main.cpp.o" \
"CMakeFiles/supra.dir/src/analyzer.cpp.o" \
"CMakeFiles/supra.dir/src/conformergenerator.cpp.o" \
"CMakeFiles/supra.dir/src/helper.cpp.o" \
"CMakeFiles/supra.dir/src/hungarian.cpp.o" \
"CMakeFiles/supra.dir/src/structure.cpp.o" \
"CMakeFiles/supra.dir/src/symmetry.cpp.o"

# External object files for target supra
supra_EXTERNAL_OBJECTS =

supra: CMakeFiles/supra.dir/main_files/main.cpp.o
supra: CMakeFiles/supra.dir/src/analyzer.cpp.o
supra: CMakeFiles/supra.dir/src/conformergenerator.cpp.o
supra: CMakeFiles/supra.dir/src/helper.cpp.o
supra: CMakeFiles/supra.dir/src/hungarian.cpp.o
supra: CMakeFiles/supra.dir/src/structure.cpp.o
supra: CMakeFiles/supra.dir/src/symmetry.cpp.o
supra: CMakeFiles/supra.dir/build.make
supra: CMakeFiles/supra.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/baum/revSUPRA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable supra"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/supra.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/supra.dir/build: supra

.PHONY : CMakeFiles/supra.dir/build

CMakeFiles/supra.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/supra.dir/cmake_clean.cmake
.PHONY : CMakeFiles/supra.dir/clean

CMakeFiles/supra.dir/depend:
	cd /home/baum/revSUPRA/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/baum/revSUPRA /home/baum/revSUPRA /home/baum/revSUPRA/build /home/baum/revSUPRA/build /home/baum/revSUPRA/build/CMakeFiles/supra.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/supra.dir/depend

