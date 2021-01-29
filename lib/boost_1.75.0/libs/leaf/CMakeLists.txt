cmake_minimum_required(VERSION 3.8)

project(leaf VERSION 0.2.1 LANGUAGES CXX)

add_library(leaf INTERFACE)

target_include_directories(leaf INTERFACE
 $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
 $<INSTALL_INTERFACE:include>
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/leaf-config-version.cmake"
  COMPATIBILITY SameMajorVersion
)

include(GNUInstallDirs)
install(TARGETS leaf
  EXPORT leaf-targets
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT leaf-targets
  FILE
    leaf-config.cmake
  NAMESPACE
    zajo::
  DESTINATION
    ${CMAKE_INSTALL_DATADIR}/cmake/leaf
)

install(FILES "${PROJECT_BINARY_DIR}/leaf-config-version.cmake"
        DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/leaf)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})