cmake_minimum_required(VERSION 3.9)
project(math)

set(CMAKE_CXX_STANDARD 11)


### googletest include
add_subdirectory(lib/gtest_1.7.0)

include_directories(
  ${gtest_SOURCE_DIR}/include
)
### googletest include

### boost include
# set (BOOST_ROOT "lib/boost_1.66.0")
# set(Boost_USE_STATIC_LIBS OFF) 
# set(Boost_USE_MULTITHREADED ON)  
# set(Boost_USE_STATIC_RUNTIME OFF) 

# FIND_PACKAGE( Boost 1.66.0 COMPONENTS program_options REQUIRED )

### boost include

add_subdirectory(lib/cvodes-3.1.0/)

include_directories(
  lib/boost_1.66.0
  lib/cvodes-3.1.0/include
  lib/eigen_3.3.3
  lib/idas-2.1.0/include
  ${CMAKE_SOURCE_DIR}
)


enable_testing()

add_executable(broken_test
        test/unit/math/rev/core/agrad_test.cpp
)

target_link_libraries (broken_test gtest gtest_main)
