cmake_minimum_required(VERSION 3.5)
project(math)

set(CMAKE_CXX_STANDARD 14)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
# Setting -O3 to work around linking issue which also affects the non-cmake build
# Dependencies aren't correctly specified, so the build only succeeds if certain
# symbols are optimized out. See #806.


if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER MATCHES "[Cc]lang")
  # This is a workaround for an old clang bug.
  # See https://bugs.llvm.org/show_bug.cgi?id=13530 for info on the clang bug.
  # The non-cmake build works around this by specifying a different libstdc++,
  # which is incorrect for many configurations, including some configurations
  # with clang < 3.9.
  if (CLANG_VERSION_STRING VERSION_LESS 3.9)
    add_compile_options("-D__STRICT_ANSI__")
  endif()
endif()

### googletest include
add_subdirectory(lib/gtest_1.7.0)

include_directories(
  ${gtest_SOURCE_DIR}/include
)

### boost include
# set (BOOST_ROOT "lib/boost_1.66.0")
# set(Boost_USE_STATIC_LIBS OFF)
# set(Boost_USE_MULTITHREADED ON)
# set(Boost_USE_STATIC_RUNTIME OFF)

# FIND_PACKAGE( Boost 1.66.0 COMPONENTS program_options REQUIRED )

### boost include

# TODO: don't build all of cvodes.
add_subdirectory(lib/cvodes-3.1.0)
# add_subdirectory(lib/idas-2.1.0)

include_directories(
  lib/boost_1.66.0
  lib/cvodes-3.1.0/include
  lib/eigen_3.3.3
  lib/idas-2.1.0/include
  ${CMAKE_SOURCE_DIR}
)

file(GLOB_RECURSE EXECUTABLES RELATIVE ${CMAKE_SOURCE_DIR} test/unit/*_test.cpp)

foreach(testsourcefile ${EXECUTABLES})
  string(REPLACE ".cpp" "" testname ${testsourcefile} )
  string(REPLACE "/" "_" testname ${testname} )
  add_executable(${testname} ${testsourcefile})
  target_link_libraries (${testname} gtest gtest_main sundials_cvodes_static
    sundials_nvecserial_static)
endforeach(testsourcefile ${EXECUTABLES})
