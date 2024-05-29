cmake_minimum_required(VERSION 3.20.2)
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
  -Wall)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(CheckIPOSupported)
check_ipo_supported(RESULT flto_support OUTPUT error LANGUAGES CXX)

if(flto_support)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(WARNING "IPO/LTO is not supported: ${error}")
endif()

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

add_subdirectory(cmake EXCLUDE_FROM_ALL)

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


add_subdirectory(test)
