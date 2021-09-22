set(TRIVIAL_SOURCES
    trivial.cpp )

set(OPENCLHPP_HEADERS
    ${CLCPP_SOURCE_DIR}/include/CL/opencl.hpp)

set(TRIVIAL_SIZET_COMPAT_SOURCES
    trivialSizeTCompat.cpp )

set(HEADEREXAMPLE_SOURCES
    headerexample.cpp )

include_directories(
    ${CLCPP_SOURCE_DIR}/include
    ${OPENCL_DIST_DIR}/include
    ${OPENCL_INCLUDE_DIR})

link_directories( ${OPENCL_LIB_DIR} )

if( MSVC )
  set(CMAKE_CXX_FLAGS " -W3 /EHsc" )
elseif( CMAKE_COMPILER_IS_GNUCXX )
  set(CMAKE_CXX_FLAGS "-Wall" )
elseif( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
  set(CMAKE_CXX_FLAGS "-Wall" )
endif()

add_executable(trivial ${TRIVIAL_SOURCES} ${OPENCLHPP_HEADERS})
target_link_libraries( trivial OpenCL )

add_executable(trivialSizeTCompat ${TRIVIAL_SIZET_COMPAT_SOURCES} ${OPENCLHPP_HEADERS})
target_link_libraries( trivialSizeTCompat OpenCL )

add_executable(headerexample ${HEADEREXAMPLE_SOURCES} ${OPENCLHPP_HEADERS})
target_link_libraries( headerexample OpenCL )
