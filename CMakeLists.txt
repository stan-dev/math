cmake_minimum_required(VERSION 3.14)
project("StanMath" VERSION 4.0.1
                 DESCRIPTION "A light-weight header-only C++ library of automatic differentiation.")

option(STANMATH_ENABLE_TEST "Enable unit tests to be built." ON)
option(STANMATH_ENABLE_COVERAGE "Build FastAD with coverage" OFF)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.10.0
)
FetchContent_MakeAvailable(googletest)

FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen
  GIT_TAG        3.3.9
)
FetchContent_MakeAvailable(eigen)


FetchContent_Declare(
  sundials
  GIT_REPOSITORY https://github.com/LLNL/sundials
  GIT_TAG        v5.7.0
)
FetchContent_MakeAvailable(sundials)


FetchContent_Declare(
  oneTBB
  GIT_REPOSITORY https://github.com/oneapi-src/oneTBB
  GIT_TAG        v2021.1.1
)
FetchContent_MakeAvailable(oneTBB)

# This is to make this library portable to other machines.
# This will be used for install.
include(GNUInstallDirs)


# Add this library as interface (header-only)
add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Set C++17 standard for project target
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_14)

# Set install destinations
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Create FastADConfigVersion.cmake which contains current project version
# This is supposed to help with (major) version compatibility.
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION ${PROJECT_VERSION}
                                 COMPATIBILITY SameMajorVersion)
configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION
    ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(EXPORT ${PROJECT_NAME}_Targets
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(FILES
        "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})

# Automate the choosing of config
# if CMAKE_BUILD_TYPE not defined
if (NOT CMAKE_BUILD_TYPE)
    # if binary directory ends with "release", use release mode
    if (${PROJECT_BINARY_DIR} MATCHES "release$")
        set(CMAKE_BUILD_TYPE RELEASE)
    # otherwise, use debug mode
    else()
        set(CMAKE_BUILD_TYPE DEBUG)
    endif()
endif()
message(STATUS "Compiling in ${CMAKE_BUILD_TYPE} mode")

if (STANMATH_ENABLE_TEST)
    # This will perform memcheck
    include(CTest)
    enable_testing()

    # Points to the root of Google Test
    set(GTEST_DIR ${PROJECT_SOURCE_DIR}/lib/googletest)

    # add test subdirectory
    add_subdirectory(${PROJECT_SOURCE_DIR}/test ${PROJECT_BINARY_DIR}/test)
endif()

# add lib subdirectory
add_subdirectory(${PROJECT_SOURCE_DIR}/lib ${PROJECT_BINARY_DIR}/lib)
