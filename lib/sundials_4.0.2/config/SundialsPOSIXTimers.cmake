# ---------------------------------------------------------------------------
# Programmer: David J. Gardner @ LLNL
# ---------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------------------
# Check if POSIX timers are available and test if they can be used as some
# compiler flags will preclude using POSIX timers even if they are present
# (e.g., -ansi).
# ---------------------------------------------------------------------------

# clock-monotonic, see if we need to link with rt
include(CheckSymbolExists)

# save and overwrite required libraries to check for timers
set(CMAKE_REQUIRED_LIBRARIES_SAVE ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES rt)

# check if _POSIX_TIMERS macro is defined in required libraries
check_symbol_exists(_POSIX_TIMERS "unistd.h;time.h" SUNDIALS_POSIX_TIMERS)

# restore required libraries
set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_SAVE})

if(SUNDIALS_POSIX_TIMERS)

  # locate rt library and hide cache variable
  find_library(SUNDIALS_RT_LIBRARY NAMES rt)
  mark_as_advanced(SUNDIALS_RT_LIBRARY)

  if(SUNDIALS_RT_LIBRARY)

    # Test timers with a simple program
    set(POSIXTest_DIR ${PROJECT_BINARY_DIR}/PosixTimersTest)
    file(MAKE_DIRECTORY ${POSIXTest_DIR})

    # Create a CMakeLists.txt file which will generate the test executable
    file(WRITE ${POSIXTest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)\n"
      "PROJECT(posixtimerstest C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
      "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(posixtimerstest posixtimerstest.c)\n"
      "TARGET_LINK_LIBRARIES(posixtimerstest ${SUNDIALS_RT_LIBRARY})\n")

    # Create a simple C source for testing
    file(WRITE ${POSIXTest_DIR}/posixtimerstest.c
      "#include <time.h>\n"
      "#include <unistd.h>\n"
      "int main(){\n"
      "time_t base_time_tv_sec = 0;\n"
      "struct timespec spec;\n"
      "clock_gettime(CLOCK_MONOTONIC_RAW, &spec);\n"
      "base_time_tv_sec = spec.tv_sec;\n"
      "clock_getres(CLOCK_MONOTONIC_RAW, &spec);\n"
      "return(0);\n"
      "}\n")

    # Use TRY_COMPILE to make the target
    try_compile(POSIX_TIMERS_TEST_OK ${POSIXTest_DIR} ${POSIXTest_DIR}
      posixtimerstest OUTPUT_VARIABLE MY_OUTPUT)

    # To ensure we do not use stuff from the previous attempts,
    # we must remove the CMakeFiles directory.
    file(REMOVE_RECURSE ${POSIXTest_DIR}/CMakeFiles)

    # Process test result
    if(POSIX_TIMERS_TEST_OK)
      # set sundials_config.h symbol
      set(SUNDIALS_HAVE_POSIX_TIMERS TRUE)
      # add rt library to list of extra libraries linked against
      set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${SUNDIALS_RT_LIBRARY})
    endif()

  endif()
endif()
