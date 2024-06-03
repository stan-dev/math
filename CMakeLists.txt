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
option(STAN_NO_RANGE_CHECKS "Disable range checks within the Stan library" OFF)
option(STAN_OPENCL "Enable OpenCL support" OFF)
option(STAN_THREADS "Enable multi-threading support" OFF)
option(STAN_MPI "Enable MPI support" OFF)

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
  -DBOOST_DISABLE_ASSERTS
  -Wall -Wno-error)

if(STAN_THREADS)
    add_compile_definitions(STAN_THREADS)
endif()

if(STAN_MPI)
    find_package(MPI REQUIRED)
    add_compile_definitions(STAN_MPI)
endif()

# Handle OpenCL if necessary
if(STAN_OPENCL)
    find_package(OpenCL REQUIRED)
    add_compile_definitions(STAN_OPENCL OPENCL_DEVICE_ID=0 OPENCL_PLATFORM_ID=0)
    # TODO: Make these per target
    include_directories(${OpenCL_INCLUDE_DIRS})
    link_libraries(${OpenCL_LIBRARIES})
endif()

# Externally provided libraries
FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG main)
FetchContent_Declare(googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG main) # need master for benchmark::benchmark

FetchContent_MakeAvailable(googletest googlebenchmark)


if(STAN_USE_SYSTEM_TBB)
  find_package(TBB REQUIRED)
else()
  FetchContent_Declare(
    tbb
    GIT_REPOSITORY https://github.com/oneapi-src/oneTBB
    GIT_TAG        v2021.7.0  # adjust this to the version you need
    CMAKE_ARGS -DTBB_STRICT=OFF -DTBB_TEST=OFF
  )
  FetchContent_MakeAvailable(tbb)
endif()
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
  FetchContent_Declare(
    boost
    URL https://boostorg.jfrog.io/artifactory/main/release/1.85.0/source/boost_1_85_0.tar.gz
  )
  FetchContent_GetProperties(boost)
  if(NOT boost_POPULATED)
    FetchContent_Populate(boost)
    set(boost_INCLUDE_DIR ${boost_SOURCE_DIR})
  endif()
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

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(CheckIPOSupported)
check_ipo_supported(RESULT flto_support OUTPUT error LANGUAGES CXX)

if(flto_support)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(WARNING "IPO/LTO is not supported: ${error}")
endif()
        
add_subdirectory(test)
