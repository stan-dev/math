# Copyright 2019, 2020 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

if(NOT CMAKE_VERSION VERSION_LESS 3.10)
  include_guard()
endif()

include(BoostMessage)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Variables

if(WIN32)
  set(__boost_default_layout "versioned")
else()
  set(__boost_default_layout "system")
endif()

set(BOOST_INSTALL_LAYOUT ${__boost_default_layout} CACHE STRING "Installation layout (versioned, tagged, or system)")
set_property(CACHE BOOST_INSTALL_LAYOUT PROPERTY STRINGS versioned tagged system)

set(BOOST_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}" CACHE STRING "Installation directory for library files")
set(BOOST_INSTALL_CMAKEDIR "${BOOST_INSTALL_LIBDIR}/cmake" CACHE STRING "Installation directory for CMake configuration files")
set(BOOST_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}" CACHE STRING "Installation directory for header files")

if(BOOST_INSTALL_LAYOUT STREQUAL "versioned")
  set(BOOST_INSTALL_INCLUDEDIR "$CACHE{BOOST_INSTALL_INCLUDEDIR}/boost-${PROJECT_VERSION_MAJOR}_${PROJECT_VERSION_MINOR}")
endif()

#

function(__boost_install_set_output_name LIB TYPE VERSION)

  set(name_debug ${LIB})
  set(name_release ${LIB})

  # prefix
  if(WIN32 AND TYPE STREQUAL "STATIC_LIBRARY")
    set_target_properties(${LIB} PROPERTIES PREFIX "lib")
  endif()

  # toolset
  if(BOOST_INSTALL_LAYOUT STREQUAL versioned)

    string(TOLOWER ${CMAKE_CXX_COMPILER_ID} toolset)

    if(toolset STREQUAL "msvc")

      set(toolset "vc")

      if(CMAKE_CXX_COMPILER_VERSION MATCHES "^([0-9]+)[.]([0-9]+)")

        if(CMAKE_MATCH_1 GREATER 18)
          math(EXPR major ${CMAKE_MATCH_1}-5)
        else()
          math(EXPR major ${CMAKE_MATCH_1}-6)
        endif()

        math(EXPR minor ${CMAKE_MATCH_2}/10)

        string(APPEND toolset ${major}${minor})

      endif()

    else()

      if(toolset STREQUAL "gnu")

        set(toolset "gcc")

      elseif(toolset STREQUAL "clang")

        if(MSVC)
          set(toolset "clangw")
        endif()

      endif()

      if(CMAKE_CXX_COMPILER_VERSION MATCHES "^([0-9]+)[.]")
        string(APPEND toolset ${CMAKE_MATCH_1})
      endif()

    endif()

    string(APPEND name_debug "-${toolset}")
    string(APPEND name_release "-${toolset}")

  endif()

  if(BOOST_INSTALL_LAYOUT STREQUAL versioned OR BOOST_INSTALL_LAYOUT STREQUAL tagged)

    # threading
    string(APPEND name_debug "-mt")
    string(APPEND name_release "-mt")

    # ABI tag

    if(MSVC)

      get_target_property(MSVC_RUNTIME_LIBRARY ${LIB} MSVC_RUNTIME_LIBRARY)

      if(MSVC_RUNTIME_LIBRARY STREQUAL "MultiThreaded$<$<CONFIG:Debug>:Debug>")

        string(APPEND name_debug "-sgd")
        string(APPEND name_release "-s")

      else()

        string(APPEND name_debug "-gd")

      endif()

    else()

      string(APPEND name_debug "-d")

    endif()

    # Arch and model
    math(EXPR bits ${CMAKE_SIZEOF_VOID_P}*8)

    string(APPEND name_debug "-x${bits}") # x86 only for now
    string(APPEND name_release "-x${bits}")

  endif()

  if(BOOST_INSTALL_LAYOUT STREQUAL versioned)

    string(REGEX REPLACE "^([0-9]+)[.]([0-9]+).*" "\\1_\\2" __ver ${VERSION})

    string(APPEND name_debug "-${__ver}")
    string(APPEND name_release "-${__ver}")

  endif()

  set_target_properties(${LIB} PROPERTIES OUTPUT_NAME_DEBUG ${name_debug})
  set_target_properties(${LIB} PROPERTIES OUTPUT_NAME ${name_release})

  if(TYPE STREQUAL "STATIC_LIBRARY")

    set_target_properties(${LIB} PROPERTIES COMPILE_PDB_NAME_DEBUG "${name_debug}")
    set_target_properties(${LIB} PROPERTIES COMPILE_PDB_NAME "${name_release}")

  endif()

endfunction()

function(__boost_install_update_include_directory lib incdir prop)

  get_target_property(value ${lib} ${prop})

  if(value STREQUAL incdir)

    set_target_properties(${lib} PROPERTIES ${prop} "$<BUILD_INTERFACE:${incdir}>;$<INSTALL_INTERFACE:${BOOST_INSTALL_INCLUDEDIR}>")

  endif()

endfunction()

# Installs a single target
# boost_install_target(TARGET target VERSION version [HEADER_DIRECTORY directory])

function(boost_install_target)

  cmake_parse_arguments(_ "" "TARGET;VERSION;HEADER_DIRECTORY" "" ${ARGN})

  if(NOT __TARGET)

    message(SEND_ERROR "boost_install_target: TARGET not given.")
    return()

  endif()

  if(NOT __VERSION)

    message(SEND_ERROR "boost_install_target: VERSION not given, but is required for installation.")
    return()

  endif()

  set(LIB ${__TARGET})

  if(NOT __HEADER_DIRECTORY)

    set(__HEADER_DIRECTORY "${PROJECT_SOURCE_DIR}/include")

  endif()

  get_target_property(TYPE ${LIB} TYPE)

  __boost_install_update_include_directory(${LIB} "${__HEADER_DIRECTORY}" INTERFACE_INCLUDE_DIRECTORIES)

  if(TYPE STREQUAL "STATIC_LIBRARY" OR TYPE STREQUAL "SHARED_LIBRARY")

    __boost_install_update_include_directory(${LIB} "${__HEADER_DIRECTORY}" INCLUDE_DIRECTORIES)

    get_target_property(OUTPUT_NAME ${LIB} OUTPUT_NAME)

    if(NOT OUTPUT_NAME)
      __boost_install_set_output_name(${LIB} ${TYPE} ${__VERSION})
    endif()

  endif()

  if(TYPE STREQUAL "SHARED_LIBRARY" OR TYPE STREQUAL "EXECUTABLE")

    get_target_property(VERSION ${LIB} VERSION)

    if(NOT VERSION)
      set_target_properties(${LIB} PROPERTIES VERSION ${__VERSION})
    endif()

  endif()

  if(LIB MATCHES "^boost_(.*)$")
    set_target_properties(${LIB} PROPERTIES EXPORT_NAME ${CMAKE_MATCH_1})
  endif()

  set(CONFIG_INSTALL_DIR "${BOOST_INSTALL_CMAKEDIR}/${LIB}-${__VERSION}")

  install(TARGETS ${LIB} EXPORT ${LIB}-targets DESTINATION ${BOOST_INSTALL_LIBDIR})

  if(WIN32 AND TYPE STREQUAL "SHARED_LIBRARY")

    install(FILES $<TARGET_PDB_FILE:${LIB}> DESTINATION ${BOOST_INSTALL_LIBDIR} OPTIONAL)

  endif()

  if(WIN32 AND TYPE STREQUAL "STATIC_LIBRARY" AND NOT CMAKE_VERSION VERSION_LESS 3.15)

    install(FILES "$<TARGET_FILE_DIR:${LIB}>/$<TARGET_FILE_PREFIX:${LIB}>$<TARGET_FILE_BASE_NAME:${LIB}>.pdb" DESTINATION ${BOOST_INSTALL_LIBDIR} OPTIONAL)

  endif()

  install(EXPORT ${LIB}-targets DESTINATION "${CONFIG_INSTALL_DIR}" NAMESPACE Boost:: FILE ${LIB}-targets.cmake)

  set(CONFIG_FILE_NAME "${CMAKE_CURRENT_BINARY_DIR}/tmpinst/${LIB}-config.cmake")
  set(CONFIG_FILE_CONTENTS "# Generated by BoostInstall.cmake for ${LIB}-${__VERSION}\n\n")

  get_target_property(INTERFACE_LINK_LIBRARIES ${LIB} INTERFACE_LINK_LIBRARIES)

  set(LINK_LIBRARIES "")

  if(TYPE STREQUAL "STATIC_LIBRARY" OR TYPE STREQUAL "SHARED_LIBRARY")
    get_target_property(LINK_LIBRARIES ${LIB} LINK_LIBRARIES)
  endif()

  if(INTERFACE_LINK_LIBRARIES OR LINK_LIBRARIES)

    string(APPEND CONFIG_FILE_CONTENTS "include(CMakeFindDependencyMacro)\n\n")

    foreach(dep IN LISTS INTERFACE_LINK_LIBRARIES LINK_LIBRARIES)

      if(${dep} MATCHES "^Boost::(.*)$")

        string(APPEND CONFIG_FILE_CONTENTS "find_dependency(boost_${CMAKE_MATCH_1} ${__VERSION} EXACT)\n")

      endif()

    endforeach()

    string(APPEND CONFIG_FILE_CONTENTS "\n")

  endif()

  string(APPEND CONFIG_FILE_CONTENTS "include(\"\${CMAKE_CURRENT_LIST_DIR}/${LIB}-targets.cmake\")\n")

  file(WRITE "${CONFIG_FILE_NAME}" "${CONFIG_FILE_CONTENTS}")
  install(FILES "${CONFIG_FILE_NAME}" DESTINATION "${CONFIG_INSTALL_DIR}")

  set(CONFIG_VERSION_FILE_NAME "${CMAKE_CURRENT_BINARY_DIR}/tmpinst/${LIB}-config-version.cmake")

  if(TYPE STREQUAL "INTERFACE_LIBRARY")

    # Header-only libraries are architecture-independent

    if(NOT CMAKE_VERSION VERSION_LESS 3.14)

      write_basic_package_version_file("${CONFIG_VERSION_FILE_NAME}" COMPATIBILITY AnyNewerVersion ARCH_INDEPENDENT)

    else()

      set(OLD_CMAKE_SIZEOF_VOID_P ${CMAKE_SIZEOF_VOID_P})
      set(CMAKE_SIZEOF_VOID_P "")

      write_basic_package_version_file("${CONFIG_VERSION_FILE_NAME}" COMPATIBILITY AnyNewerVersion)

      set(CMAKE_SIZEOF_VOID_P ${OLD_CMAKE_SIZEOF_VOID_P})

    endif()

  else()

    write_basic_package_version_file("${CONFIG_VERSION_FILE_NAME}" COMPATIBILITY AnyNewerVersion)

  endif()

  install(FILES "${CONFIG_VERSION_FILE_NAME}" DESTINATION "${CONFIG_INSTALL_DIR}")

endfunction()

# boost_install([VERSION version] [TARGETS targets...] [HEADER_DIRECTORY directory])

function(boost_install)

  cmake_parse_arguments(_ "" "VERSION;HEADER_DIRECTORY" "TARGETS" ${ARGN})

  if(NOT __VERSION)

    if(NOT PROJECT_VERSION)

      message(AUTHOR_WARNING "boost_install: VERSION not given, PROJECT_VERSION not set, but a version is required for installation.")
      return()

    else()

      boost_message(DEBUG "boost_install: VERSION not given, assuming PROJECT_VERSION ('${PROJECT_VERSION}')")
      set(__VERSION ${PROJECT_VERSION})

    endif()

  endif()

  if(__UNPARSED_ARGUMENTS)

    message(AUTHOR_WARNING "boost_install: extra arguments ignored: ${__UNPARSED_ARGUMENTS}")

  endif()

  if(__HEADER_DIRECTORY)

    get_filename_component(__HEADER_DIRECTORY "${__HEADER_DIRECTORY}" ABSOLUTE)
    install(DIRECTORY "${__HEADER_DIRECTORY}/" DESTINATION "${BOOST_INSTALL_INCLUDEDIR}")

  endif()

  foreach(target IN LISTS __TARGETS)

    boost_install_target(TARGET ${target} VERSION ${__VERSION} HEADER_DIRECTORY ${__HEADER_DIRECTORY})

  endforeach()

endfunction()
