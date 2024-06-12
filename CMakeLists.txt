cmake_minimum_required(VERSION 3.16.3)
project(
    stanmath
    VERSION 0.0.1
    LANGUAGES C CXX)

include(FetchContent)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)
set(CMAKE_VERBOSE_MAKEFILE YES)
cmake_policy(SET CMP0069 NEW)
# Configuration Options
option(STAN_BUILD_DOCS "Build the Stan Math library documentation" OFF)
option(STAN_TEST_HEADERS "Build the targets for the header checks" OFF)
option(STAN_NO_RANGE_CHECKS "Disable range checks within the Stan library" OFF)
option(STAN_MPI "Enable MPI support" OFF)
option(STAN_OPENCL "Enable OpenCL support" OFF)
set(STAN_OPENCL_DEVICE_ID "0" CACHE STRING "Set the OpenCL Device ID at compile time" FORCE)
set(STAN_OPENCL_PLATFORM_ID "0" CACHE STRING "Set the OpenCL Platform ID at compile time" FORCE)
set(STAN_INTEGRATED_OPENCL "0" CACHE STRING "Whether the platform and device can use integrated opencl features" FORCE)
message(STATUS "OpenCL Platform: " ${STAN_OPENCL_PLATFORM_ID})
message(STATUS "OpenCL Device: " ${STAN_OPENCL_DEVICE_ID})
option(STAN_THREADS "Enable multi-threading support" OFF)
if(STAN_NO_RANGE_CHECKS)
    add_compile_definitions(STAN_NO_RANGE_CHECKS)
endif()
if(POLICY CMP0069)
  cmake_policy(SET CMP0069 NEW)
endif()

# Set compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wno-deprecated-declarations)
    if(APPLE)
        add_compile_options(-Wno-unknown-warning-option -Wno-tautological-compare -Wno-sign-compare)
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wno-sign-compare)
endif()

add_compile_options(
  -DNO_FPRINTF_OUTPUT
  -DBOOST_DISABLE_ASSERTS
  -DTBB_INTERFACE_NEW
  -D_REENTRANT
  -Wno-deprecated-declarations
  -Wall )

if(STAN_THREADS)
    add_compile_definitions(STAN_THREADS)
endif()

if(STAN_MPI)
    find_package(MPI REQUIRED)
    add_compile_definitions(STAN_MPI)
    add_compile_options(-Wno-delete-non-virtual-dtor)
endif()

# Handle OpenCL if necessary
if(STAN_OPENCL)
    # Externally provided libraries
    FetchContent_Declare(OpenCLHeaders
            GIT_REPOSITORY https://github.com/KhronosGroup/OpenCL-CLHPP
            GIT_TAG v2.0.15)
    FetchContent_MakeAvailable(OpenCLHeaders)
    find_package(OpenCL REQUIRED)
    add_compile_definitions(STAN_OPENCL OPENCL_DEVICE_ID=${STAN_OPENCL_DEVICE_ID}
      OPENCL_PLATFORM_ID=${STAN_OPENCL_PLATFORM_ID} CL_HPP_TARGET_OPENCL_VERSION=120
      CL_HPP_MINIMUM_OPENCL_VERSION=120 CL_HPP_ENABLE_EXCEPTIONS INTEGRATED_OPENCL=${STAN_INTEGRATED_OPENCL})
    add_compile_options(-Wno-ignored-attributes)

endif()

# Externally provided libraries
FetchContent_Declare(googletest
        DOWNLOAD_EXTRACT_TIMESTAMP ON
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG main)
FetchContent_Declare(googlebenchmark
        DOWNLOAD_EXTRACT_TIMESTAMP ON
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG main) # need master for benchmark::benchmark

FetchContent_MakeAvailable(googletest googlebenchmark)


if(STAN_USE_SYSTEM_TBB)
  find_package(TBB REQUIRED)
else()
  FetchContent_Declare(
    tbb
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    GIT_REPOSITORY https://github.com/oneapi-src/oneTBB
    GIT_TAG        v2021.7.0  # adjust this to the version you need
    CMAKE_ARGS -DTBB_STRICT=OFF -DTBB_TEST=OFF
    CMAKE_CACHE_ARGS -DTBB_STRICT:BOOL=OFF -DTBB_TEST:BOOL=OFF
  )
  FetchContent_GetProperties(tbb)
    if(NOT tbb_POPULATED)
    FetchContent_Populate(tbb)
    # Add the fetched content to your build
    set(TBB_STRICT OFF CACHE BOOL "Treat compiler warnings as errors")
    add_subdirectory(${tbb_SOURCE_DIR} ${tbb_BINARY_DIR})
    endif()
endif()
# For tbb
add_compile_options(-Wno-error -Wno-unused-value)
if(POLICY CMP0069)
  cmake_policy(SET CMP0069 NEW)
endif()
if(STAN_USE_SYSTEM_EIGEN)
  find_package(Eigen3 REQUIRED)
else()
  set(EIGEN_BUILD_DOC OFF)
  # note: To disable eigen tests,
  # you should put this code in a add_subdirectory to avoid to change
  # BUILD_TESTING for your own project too since variables are directory
  # scoped
  set(EIGEN_BUILD_BTL OFF)
  set(EIGEN_BUILD_TESTING OFF)
  set(EIGEN_BUILD_PKGCONFIG OFF)
  set(EIGEN_LEAVE_TEST_IN_ALL_TARGET OFF)
  FetchContent_Declare(
    Eigen
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0)
  FetchContent_MakeAvailable(Eigen)
  #TODO: Need to just put dependency file in subdirectory
  # since macros are directory scoped
  set(BUILD_TESTING ON)
  set(EIGEN_BUILD_TESTING ON)
endif()

if (STAN_USE_SYSTEM_SUNDIALS)
  find_package(SUNDIALS REQUIRED)
else()
  FetchContent_Declare(
    sundials
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    GIT_REPOSITORY https://github.com/LLNL/sundials
    GIT_TAG        v6.1.1
      # adjust this to the version you need
  )
  FetchContent_GetProperties(sundials)
  if(NOT sundials_POPULATED)
    FetchContent_Populate(sundials)
    add_subdirectory(${sundials_SOURCE_DIR} ${sundials_BINARY_DIR})
  endif()
endif()


if (STAN_USE_SYSTEM_BOOST)
  find_package(Boost REQUIRED)
else()
  set(Boost_USE_STATIC_LIBS ON)
  set(Boost_USE_STATIC_RUNTIME ON)
  set(BOOST_ENABLE_CMAKE ON)
  set(BUILD_SHARED_LIBS OFF)
  set(BOOST_DETAILED_CONFIGURE ON)
  set(BOOST_ENABLE_PYTHON OFF)
  set(BOOST_BUILD_TESTS OFF)
  set(CMAKE_BUILD_TYPE Release)

  set(BOOST_INCLUDE_LIBRARIES math numeric/odeint lexical_cast optional random)
  if (STAN_MPI)
    set(BOOST_ENABLE_MPI ON)
    set(BOOST_MPI ON)
    set(BOOST_INCLUDE_LIBRARIES ${BOOST_INCLUDE_LIBRARIES} mpi serialization)
  endif()
  FetchContent_Declare(
    Boost
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    URL https://github.com/boostorg/boost/releases/download/boost-1.85.0/boost-1.85.0-cmake.tar.xz
  )
  FetchContent_MakeAvailable(Boost)
endif()
# Library target
add_library(stanmath INTERFACE)  # Assuming you only have headers
target_include_directories(stanmath INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)

# If you have sources, specify them and use add_library(stanmath SHARED or STATIC) instead

include(GNUInstallDirs)
install(DIRECTORY stan/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/stan)
install(TARGETS stanmath
        EXPORT stanmathTargets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Export targets
install(EXPORT stanmathTargets
        FILE stanmathTargets.cmake
        NAMESPACE stanmath::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/stanmath)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "stanmathConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

install(FILES "stanmathConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/stanmathConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/stanmath)

# add docs
if (STAN_BUILD_DOCS)
find_package(Doxygen REQUIRED)

add_custom_target(
  doxygen ALL
  COMMAND mkdir -p doc/api
  COMMAND
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/doxygen.cfg
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/pretty_stuff/eigen_navtree_hacks.js ${CMAKE_CURRENT_SOURCE_DIR}/doc/api/html/eigen_navtree_hacks.js
  COMMENT "Generating API documentation with Doxygen"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(CheckIPOSupported)
check_ipo_supported(RESULT flto_support OUTPUT error LANGUAGES CXX)

# Becuase of bug with mingw
# lto1.exe: error: two or more sections for .gnu.lto__ZN4stan...
IF (NOT WIN32)
    if(flto_support)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(WARNING "IPO/LTO is not supported: ${error}")
    endif()
endif()
add_subdirectory(test)


