cmake_minimum_required(VERSION 3.5)
project(math)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3")


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
