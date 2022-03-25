# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# Test if POSIX timers are available and can be used (some compiler flags will
# disable POSIX timers even if they are present e.g., -ansi).
#
# Note -D_POSIX_C_SOURCE=199309L is the minimum POSIX version needed for struct
# timespec and clock_monotonic()
# ------------------------------------------------------------------------------

macro(posix_timers_test)

  set(options )
  set(oneValueArgs POSIX RT_LIB)
  set(multiValueArgs )

  # parse keyword arguments/options
  cmake_parse_arguments(posix_timers_test
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Test timers with a simple program
  set(POSIX_TIMER_TEST_DIR ${PROJECT_BINARY_DIR}/POSIX_TIMER_TEST)
  file(MAKE_DIRECTORY ${POSIX_TIMER_TEST_DIR})

  # Create a CMakeLists.txt file which will generate the test executable
  file(WRITE ${POSIX_TIMER_TEST_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
    "PROJECT(ltest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
    "SET(CMAKE_C_STANDARD ${CMAKE_C_STANDARD})\n"
    "SET(CMAKE_C_EXTENSIONS ${CMAKE_C_EXTENSIONS})\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_COMPILE_DEFINITIONS(ltest PRIVATE \"${posix_timers_test_POSIX}\")\n"
    "TARGET_LINK_LIBRARIES(ltest \"${posix_timers_test_RT_LIB}\")\n")

  # Create a simple C source for testing
  file(WRITE ${POSIX_TIMER_TEST_DIR}/ltest.c
    "#include <time.h>\n"
    "#include <unistd.h>\n"
    "int main(){\n"
    "struct timespec spec;\n"
    "clock_gettime(CLOCK_MONOTONIC, &spec);\n"
    "clock_getres(CLOCK_MONOTONIC, &spec);\n"
    "return(0);\n"
    "}\n")

  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${POSIX_TIMER_TEST_DIR}/CMakeFiles)

  # Use TRY_COMPILE to make the target
  try_compile(COMPILE_OK ${POSIX_TIMER_TEST_DIR} ${POSIX_TIMER_TEST_DIR} ltest
    OUTPUT_VARIABLE COMPILE_OUTPUT)

endmacro()


if (NOT SUNDIALS_POSIX_TIMERS)

  # Test for timers without any modifications
  posix_timers_test()
  if(COMPILE_OK)
    message(STATUS "Looking for POSIX timers... found")
  endif()

  # Test failed, try again with -D_POSIX_C_SOURCE=199309L
  if(NOT COMPILE_OK)
    posix_timers_test(POSIX "_POSIX_C_SOURCE=199309L")
    if(COMPILE_OK)
      message(STATUS "Looking for POSIX timers (setting _POSIX_C_SOURCE)... found")
      set(POSIX_TIMERS_NEED_POSIX_C_SOURCE TRUE)
    endif()
  endif()

  # Test failed, try again linking to rt
  if(NOT COMPILE_OK)

    # Locate rt library and hide cache variable
    find_library(SUNDIALS_RT_LIBRARY NAMES rt)
    mark_as_advanced(SUNDIALS_RT_LIBRARY)

    if(SUNDIALS_RT_LIBRARY)
      message(STATUS "Looking for rt library... found")
      posix_timers_test(RT_LIB "${SUNDIALS_RT_LIBRARY}")
      if(COMPILE_OK)
        message(STATUS "Looking for POSIX timers (linking to rt)... found")
        set(POSIX_TIMERS_NEED_RT_LIBRARY TRUE)
        set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${SUNDIALS_RT_LIBRARY})
      endif()
    else()
      message(STATUS "Looking for rt library... FAILED")
    endif()

  endif()

  # Test failed, try again linking to rt and with -D_POSIX_C_SOURCE=199309L
  if((NOT COMPILE_OK) AND SUNDIALS_RT_LIBRARY)
    posix_timers_test(POSIX "_POSIX_C_SOURCE=199309L" RT_LIB "${SUNDIALS_RT_LIBRARY}")
    if(COMPILE_OK)
      message(STATUS "Looking for POSIX timers (setting _POSIX_C_SOURCE and linking to rt)... found")
      set(POSIX_TIMERS_NEED_POSIX_C_SOURCE TRUE)
      set(POSIX_TIMERS_NEED_RT_LIBRARY TRUE)
      set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${SUNDIALS_RT_LIBRARY})
    endif()
  endif()

  # Set POSIX timer status
  if(COMPILE_OK)
    set(SUNDIALS_POSIX_TIMERS TRUE)
  else()
    message(STATUS "Looking for POSIX timers... FAILED")
    set(SUNDIALS_POSIX_TIMERS FALSE)
  endif()

endif()
