add_custom_command(
    OUTPUT stripped/cl.h
    COMMAND ${CMAKE_COMMAND} -E make_directory stripped
    COMMAND ${CMAKE_COMMAND} -E make_directory mocks
    COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/strip_defines.py ${OPENCL_INCLUDE_DIR}/CL/cl.h stripped/cl.h
    COMMENT "Stripping defines from cl.h"
    DEPENDS ${OPENCL_INCLUDE_DIR}/CL/cl.h strip_defines.py)
	
add_custom_command(
    OUTPUT mocks/Mockcl.c mocks/Mockcl.h
    COMMAND ruby ${CMOCK_DIR}/lib/cmock.rb -o${CMAKE_CURRENT_SOURCE_DIR}/cmock.yml stripped/cl.h
    COMMENT "Generating mocks"
    DEPENDS stripped/cl.h cmock.yml)
	
add_custom_command(
    OUTPUT test_openclhpp_Runner.c
    COMMAND ruby ${UNITY_DIR}/auto/generate_test_runner.rb test_openclhpp.cpp cmock.yml ${CMAKE_CURRENT_BINARY_DIR}/test_openclhpp_Runner.c
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating test runner"
    DEPENDS test_openclhpp.cpp cmock.yml)
	
include_directories(
    ${CLCPP_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/mocks
    ${OPENCL_DIST_DIR}/include
    ${OPENCL_INCLUDE_DIR}
    ${OPENCL_INCLUDE_DIR}/CL
    ${UNITY_DIR}/src
    ${CMOCK_DIR}/src)

if(CMAKE_COMPILER_IS_GNUCC OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "(Apple)?Clang")
    add_compile_options(-Wno-deprecated-declarations)
endif()

add_definitions(-DCL_TARGET_OPENCL_VERSION=300)
add_definitions(-DCL_EXPERIMENTAL)

add_definitions(-DUNITY_SUPPORT_64)
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    add_definitions(-DUNITY_POINTER_WIDTH=64)
    add_definitions("-DCMOCK_MEM_PTR_AS_INT=unsigned long long")
    add_definitions(-DCMOCK_MEM_ALIGN=3)
endif()
if( CMAKE_SIZEOF_LONG EQUAL 8 )
    add_definitions(-DUNITY_LONG_WIDTH=64)
endif()

set(TEST_HEADERS
    ${CLCPP_SOURCE_DIR}/include/CL/opencl.hpp
    mocks/Mockcl.h)

set(TEST_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/test_openclhpp_Runner.c
    test_openclhpp.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/mocks/Mockcl.c
    ${CMOCK_DIR}/src/cmock.c
    ${UNITY_DIR}/src/unity.c)

# TODO enable testing for OpenCL 1.0 and 1.1
foreach(VERSION 120 200 210 220 300)
  set(TEST_EXE test_openclhpp_${VERSION})
  add_executable(${TEST_EXE} ${TEST_SOURCES} ${TEST_HEADERS})
  target_compile_definitions(${TEST_EXE}
	  PUBLIC -DCL_HPP_TARGET_OPENCL_VERSION=${VERSION}
  )
  add_test(NAME ${TEST_EXE} COMMAND ${TEST_EXE})
endforeach(VERSION)

