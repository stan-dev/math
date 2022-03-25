# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Find module that locates the MAGMA linear algebra library.
# -----------------------------------------------------------------------------

# find the MAGMA include path
find_path(MAGMA_INCLUDE_DIR magma_v2.h
  NAMES magma_v2.h
  HINTS ${MAGMA_DIR} $ENV{MAGMA_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with MAGMA header"
)

# find the main MAGMA library
find_library(MAGMA_LIBRARY
  NAMES magma
  HINTS ${MAGMA_DIR} $ENV{MAGMA_DIR}
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
  DOC "The MAGMA library.")

# Find the optional sparse component
if("SPARSE" IN_LIST MAGMA_FIND_COMPONENTS)
  set(_sparse_required MAGMA_SPARSE_LIBRARY)
  find_library(MAGMA_SPARSE_LIBRARY
    NAMES magma_sparse
    HINTS ${MAGMA_DIR} $ENV{MAGMA_DIR}
    PATH_SUFFIXES lib lib64
    NO_DEFAULT_PATH
    DOC "The MAGMA sparse library.")
else()
  set(_sparse_required )
endif()

# Determine MAGMA version and libraries it depends on
if(MAGMA_LIBRARY AND MAGMA_INCLUDE_DIR)

  get_filename_component(libdir ${MAGMA_LIBRARY} DIRECTORY)
  find_file(MAGMA_PKG_CONFIG_PATH magma.pc PATHS "${libdir}/pkgconfig")

  if(MAGMA_PKG_CONFIG_PATH)

    file(STRINGS ${MAGMA_PKG_CONFIG_PATH} _version_string REGEX "Version: [0-9].[0-9].[0-9]")
    string(REGEX MATCHALL "[0-9]" _version_full "${_version_string}")

    list(GET _version_full 0 _version_major)
    list(GET _version_full 1 _version_minor)
    list(GET _version_full 2 _version_patch)

    set(MAGMA_VERSION "${_version_major}.${_version_minor}.${_version_patch}")

    file(STRINGS ${MAGMA_PKG_CONFIG_PATH} _libraries_string REGEX "Libs:.*")
    string(REPLACE " " ";" _libraries_list ${_libraries_string})
    list(SUBLIST _libraries_list 1 -1 _libraries_list) # remove 'Libs:' part

    set(_interface_libraires )
    foreach(lib ${_libraries_list})
      if(NOT (lib STREQUAL "-lmagma" OR lib STREQUAL "-lmagma_sparse"
            OR lib STREQUAL "-L\${libdir}" OR lib STREQUAL "") )

        # Remove -l only from the beginning of the string
        string(REPLACE "^-l" "" lib ${lib})
        list(APPEND _interface_libraires ${lib})

        # Check if we need to find roc::hipblas or roc::hipsparse
        if(SUNDIALS_MAGMA_BACKENDS MATCHES "HIP")
          if((lib STREQUAL "roc::hipblas") AND (NOT TARGET roc::hipblas))
            find_package(hipblas REQUIRED)
          endif()
          if((lib STREQUAL "roc::hipsparse") AND (NOT TARGET roc::hipsparse))
            find_package(hipsparse REQUIRED)
          endif()
        endif()

      endif()
    endforeach()

  endif()
endif()

set(MAGMA_LIBRARIES "${MAGMA_LIBRARY};${_interface_libraires}")

find_package_handle_standard_args(MAGMA
  REQUIRED_VARS
    MAGMA_LIBRARY
    MAGMA_LIBRARIES
    MAGMA_INCLUDE_DIR
    ${_sparse_required}
  VERSION_VAR
    MAGMA_VERSION
  )

# Create target for MAGMA
if(MAGMA_FOUND)

  if(NOT TARGET SUNDIALS::MAGMA)
    add_library(SUNDIALS::MAGMA UNKNOWN IMPORTED)
  endif()

  set_target_properties(SUNDIALS::MAGMA PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${MAGMA_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${_interface_libraires}"
    IMPORTED_LOCATION "${MAGMA_LIBRARY}")

  if(MAGMA_SPARSE_LIBRARY)
    if(NOT TARGET SUNDIALS::MAGMA_SPARSE)
      add_library(SUNDIALS::MAGMA_SPARSE UNKNOWN IMPORTED)
    endif()

    set_target_properties(SUNDIALS::MAGMA_SPARSE PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${MAGMA_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES "${MAGMA_LIBRARY};${_interface_libraires}"
      IMPORTED_LOCATION "${MAGMA_SPARSE_LIBRARY}")
  endif()

endif()
