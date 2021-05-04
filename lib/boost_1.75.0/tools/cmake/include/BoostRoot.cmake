# Copyright 2019, 2020 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

if(CMAKE_SOURCE_DIR STREQUAL Boost_SOURCE_DIR AND WIN32 AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

    set(CMAKE_INSTALL_PREFIX "C:/Boost" CACHE PATH "Installation path prefix, prepended to installation directories" FORCE)

endif()

if(NOT BOOST_ENABLE_CMAKE)

  message(FATAL_ERROR
    "CMake support in Boost is experimental and part of an ongoing "
    "development effort. It's not ready for use yet. Please use b2 "
    "(Boost.Build) to build and install Boost.")

endif()

include(BoostMessage)
include(BoostInstall)

# --with-<library>
set(BOOST_INCLUDE_LIBRARIES "" CACHE STRING "List of libraries to build (default: all but excluded and incompatible)")

# --without-<library>
set(BOOST_EXCLUDE_LIBRARIES "" CACHE STRING "List of libraries to exclude from build")

set(BOOST_INCOMPATIBLE_LIBRARIES beast;callable_traits;compute;gil;hana;hof;safe_numerics;serialization;static_string;stl_interfaces;yap CACHE STRING "List of libraries with incompatible CMakeLists.txt files")

# --layout, --libdir, --cmakedir, --includedir in BoostInstall

# runtime-link=static|shared

set(BOOST_RUNTIME_LINK shared CACHE STRING "Runtime library selection for the MS ABI (shared or static)")
set_property(CACHE BOOST_RUNTIME_LINK PROPERTY STRINGS shared static)

if(NOT CMAKE_MSVC_RUNTIME_LIBRARY)

  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

  if(NOT BOOST_RUNTIME_LINK STREQUAL "static")
    string(APPEND CMAKE_MSVC_RUNTIME_LIBRARY "DLL")
  endif()

endif()

# Functions

function(__boost_auto_install __boost_lib)
  if(NOT CMAKE_VERSION VERSION_LESS 3.13)

    string(MAKE_C_IDENTIFIER "${__boost_lib}" __boost_lib_target)

    if(TARGET "Boost::${__boost_lib_target}" AND TARGET "boost_${__boost_lib_target}")

      get_target_property(__boost_lib_incdir "boost_${__boost_lib_target}" INTERFACE_INCLUDE_DIRECTORIES)

      if(__boost_lib_incdir STREQUAL "${BOOST_SUPERPROJECT_SOURCE_DIR}/libs/${__boost_lib}/include")

        boost_message(DEBUG "Enabling installation for '${__boost_lib}'")
        boost_install(TARGETS "boost_${__boost_lib_target}" VERSION "${BOOST_SUPERPROJECT_VERSION}" HEADER_DIRECTORY "${BOOST_SUPERPROJECT_SOURCE_DIR}/libs/${__boost_lib}/include")

      else()
        boost_message(DEBUG "Not enabling installation for '${__boost_lib}'; interface include directory '${__boost_lib_incdir}' does not equal '${BOOST_SUPERPROJECT_SOURCE_DIR}/libs/${__boost_lib}/include'")
      endif()

    else()
      boost_message(DEBUG "Not enabling installation for '${__boost_lib}'; targets 'Boost::${__boost_lib_target}' and 'boost_${__boost_lib_target}' weren't found")
    endif()

  endif()
endfunction()

function(__boost_scan_dependencies lib var)

  file(STRINGS "${BOOST_SUPERPROJECT_SOURCE_DIR}/libs/${lib}/CMakeLists.txt" data)

  set(result "")

  foreach(line IN LISTS data)

    if(line MATCHES "^[ ]*Boost::([A-Za-z0-9_]+)[ ]*$")

      string(REGEX REPLACE "^numeric_" "numeric/" dep ${CMAKE_MATCH_1})
      list(APPEND result ${dep})

    endif()

  endforeach()

  set(${var} ${result} PARENT_SCOPE)

endfunction()

#

if(CMAKE_SOURCE_DIR STREQUAL Boost_SOURCE_DIR)

  include(CTest)
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)

  # link=static|shared
  option(BUILD_SHARED_LIBS "Build shared libraries")

  # --stagedir
  set(BOOST_STAGEDIR "${CMAKE_CURRENT_BINARY_DIR}/stage" CACHE STRING "Build output directory")

  if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${BOOST_STAGEDIR}/bin")
  endif()

  if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${BOOST_STAGEDIR}/lib")
  endif()

  if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${BOOST_STAGEDIR}/lib")
  endif()

endif()

file(GLOB __boost_libraries RELATIVE "${BOOST_SUPERPROJECT_SOURCE_DIR}/libs" "${BOOST_SUPERPROJECT_SOURCE_DIR}/libs/*/CMakeLists.txt" "${BOOST_SUPERPROJECT_SOURCE_DIR}/libs/numeric/*/CMakeLists.txt")

# Check for mistakes in BOOST_INCLUDE_LIBRARIES

foreach(__boost_included_lib IN LISTS BOOST_INCLUDE_LIBRARIES)

  if(NOT "${__boost_included_lib}/CMakeLists.txt" IN_LIST __boost_libraries)

    message(WARNING "Library '${__boost_included_lib}' given in BOOST_INCLUDE_LIBRARIES has not been found.")

  endif()

endforeach()

# Scan for dependencies

set(__boost_include_libraries ${BOOST_INCLUDE_LIBRARIES})

if(__boost_include_libraries)
  list(REMOVE_DUPLICATES __boost_include_libraries)
endif()

set(__boost_libs_to_scan ${__boost_include_libraries})

while(__boost_libs_to_scan)

  boost_message(DEBUG "Scanning dependencies: ${__boost_libs_to_scan}")

  set(__boost_dependencies "")

  foreach(__boost_lib IN LISTS __boost_libs_to_scan)

    __boost_scan_dependencies(${__boost_lib} __boost_deps)
    list(APPEND __boost_dependencies ${__boost_deps})

  endforeach()

  list(REMOVE_DUPLICATES __boost_dependencies)

  set(__boost_libs_to_scan ${__boost_dependencies})

  if(__boost_libs_to_scan)
    list(REMOVE_ITEM __boost_libs_to_scan ${__boost_include_libraries})
  endif()

  list(APPEND __boost_include_libraries ${__boost_libs_to_scan})

endwhile()

# Installing targets created in other directories requires CMake 3.13
if(CMAKE_VERSION VERSION_LESS 3.13)

  boost_message(VERBOSE "Boost installation support is limited on CMake ${CMAKE_VERSION} (need 3.13)")

endif()

foreach(__boost_lib_cml IN LISTS __boost_libraries)

  get_filename_component(__boost_lib "${__boost_lib_cml}" DIRECTORY)

  if(__boost_lib IN_LIST BOOST_INCOMPATIBLE_LIBRARIES)

    boost_message(DEBUG "Skipping incompatible Boost library ${__boost_lib}")

  elseif(__boost_lib IN_LIST BOOST_EXCLUDE_LIBRARIES)

    boost_message(DEBUG "Skipping excluded Boost library ${__boost_lib}")

  elseif(NOT BOOST_INCLUDE_LIBRARIES OR __boost_lib IN_LIST BOOST_INCLUDE_LIBRARIES)

    boost_message(VERBOSE "Adding Boost library ${__boost_lib}")
    add_subdirectory(libs/${__boost_lib})

    __boost_auto_install(${__boost_lib})

  elseif(__boost_lib IN_LIST __boost_include_libraries)

    set(BUILD_TESTING OFF) # hide cache variable

    boost_message(VERBOSE "Adding dependent Boost library ${__boost_lib}")
    add_subdirectory(libs/${__boost_lib})

    __boost_auto_install(${__boost_lib})

    unset(BUILD_TESTING)

  else()

    set(BUILD_TESTING OFF) # hide cache variable

    boost_message(DEBUG "Adding Boost library ${__boost_lib} with EXCLUDE_FROM_ALL")
    add_subdirectory(libs/${__boost_lib} EXCLUDE_FROM_ALL)

    unset(BUILD_TESTING)

  endif()

endforeach()
