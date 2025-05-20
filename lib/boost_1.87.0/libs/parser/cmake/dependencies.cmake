# Copyright Louis Dionne 2016
# Copyright Zach Laine 2024
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)

###############################################################################
# Boost
###############################################################################
set(Boost_USE_STATIC_LIBS ON)

if (NOT BOOST_BRANCH)
  set(BOOST_BRANCH master)
endif()

add_custom_target(boost_clone_superproject
  DEPENDS
    ${CMAKE_BINARY_DIR}/boost_root/LICENSE_1_0.txt
  COMMENT
    "Cloning Boost superproject."
  VERBATIM)

add_custom_command(
  OUTPUT ${CMAKE_BINARY_DIR}/boost_root/LICENSE_1_0.txt
  COMMAND git clone --depth 100 -b ${BOOST_BRANCH}
    https://github.com/boostorg/boost.git boost_root
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# MSVC determines the compiler and not the target platform
# here we need to check WIN32
if (WIN32)
  set(b2_exe b2.exe)
else()
  set(b2_exe b2)
endif()

add_custom_target(boost_clone_deps
  DEPENDS
    ${CMAKE_BINARY_DIR}/boost_root/${b2_exe}
  COMMENT
    "Cloning Boost dependencies."
  VERBATIM)
add_dependencies(boost_clone_deps boost_clone_superproject)


if (WIN32)
  # Windows:
  # in order to build 'b2_exe' with the same toolset as configured in the current cmake run,
  # we need to tell the boost bootstrap process and pass some parameters
  # - as first parameter msvc, gcc or clang for the compiler type
  # - additionally, for non-MSVC we need to set the environment variable CXX to the compiler: ${CMAKE_CXX_COMPILER} which
  #   is a fully qualified path name, so that the compiler can be found by the bootstrap_cmd
  # calling b2.exe should work also for non-MSVC compilers, as during the build, the runtime directory is part of the search PATH

  if (MSVC OR CMAKE_CXX_SIMULATE_ID STREQUAL "MSVC")
    message(STATUS "configure BOOST for Visual Studio built-in compilers (i.e cl, clang-cl and clang")
    set(bootstrap_cmd ./bootstrap.bat msvc)
    # here we do not need to distinguish the different compilers as only the frontend is different
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(STATUS "configure BOOST for Clang compiler")
    set(bootstrap_cmd ./bootstrap.bat clang)
    set(COMMAND_ENV set CXX=${CMAKE_CXX_COMPILER})
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "configure BOOST for GCC compiler")
    set(bootstrap_cmd ./bootstrap.bat gcc)
    set(COMMAND_ENV set CXX=${CMAKE_CXX_COMPILER})
  endif ()
else()
  # windres produces relocations that are rejected
  # by stricter ld configurations used in some distros
  set(bootstrap_cmd env B2_DONT_EMBED_MANIFEST=true ./bootstrap.sh)
endif()

add_custom_command(
  OUTPUT ${CMAKE_BINARY_DIR}/boost_root/${b2_exe}
  COMMAND git submodule init libs/assert
  COMMAND git submodule init libs/config
  COMMAND git submodule init libs/core
  COMMAND git submodule init libs/hana
  COMMAND git submodule init tools/build
  COMMAND git submodule init libs/headers
  COMMAND git submodule init tools/boost_install
  COMMAND git submodule update --jobs 3 --depth 100
  COMMAND ${COMMAND_ENV}
  COMMAND ${bootstrap_cmd}
  COMMAND ./b2 headers
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/boost_root)

add_library(boost INTERFACE)
add_dependencies(boost boost_clone_deps)
target_include_directories(boost INTERFACE ${CMAKE_BINARY_DIR}/boost_root)
