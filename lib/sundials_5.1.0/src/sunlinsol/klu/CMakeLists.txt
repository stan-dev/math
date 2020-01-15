# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the klu SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_KLU\n\")")

# Source files for the library
set(sunlinsolklu_SOURCES sunlinsol_klu.c)

# Common SUNDIALS sources included in the library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c)

# Exported header files
set(sunlinsolklu_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_klu.h)

# Rules for building and installing the static library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Install the library
if(BUILD_STATIC_LIBS)

  add_library(sundials_sunlinsolklu_static
    STATIC ${sunlinsolklu_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolklu_static
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsolklu
    CLEAN_DIRECT_OUTPUT 1)

  # depends on sunmatrixsparse and KLU
  target_link_libraries(sundials_sunlinsolklu_static
    PUBLIC sundials_sunmatrixsparse_static ${KLU_LIBRARIES})

  target_compile_definitions(sundials_sunlinsolklu_static
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolklu_static
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the library
if(BUILD_SHARED_LIBS)

  add_library(sundials_sunlinsolklu_shared
    SHARED ${sunlinsolklu_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsolklu_shared
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsolklu
    CLEAN_DIRECT_OUTPUT 1
    VERSION ${sunlinsollib_VERSION}
    SOVERSION ${sunlinsollib_SOVERSION})

  # depends on sunmatrixsparse and KLU
  target_link_libraries(sundials_sunlinsolklu_shared
    PUBLIC sundials_sunmatrixsparse_shared ${KLU_LIBRARIES})

  target_compile_definitions(sundials_sunlinsolklu_shared
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsolklu_shared
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_SHARED_LIBS)

# Install the header files
install(FILES ${sunlinsolklu_HEADERS} DESTINATION include/sunlinsol)

#
message(STATUS "Added SUNLINSOL_KLU module")

# If FCMIX is enabled, build and install the Fortran77 library
if(F77_INTERFACE_ENABLE AND F77_FOUND)

  set(fsunlinsolklu_SOURCES fsunlinsol_klu.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fsunlinsolklu_static
      STATIC ${fsunlinsolklu_SOURCES})
    set_target_properties(sundials_fsunlinsolklu_static
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsolklu
      CLEAN_DIRECT_OUTPUT 1)
    # depends on fnvecserial, fsunmatrixsparse, sunlinsolklu
    target_link_libraries(sundials_fsunlinsolklu_static
      PUBLIC
      sundials_fnvecserial_static
      sundials_fsunmatrixsparse_static
      sundials_sunlinsolklu_static)
    install(TARGETS sundials_fsunlinsolklu_static
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fsunlinsolklu_shared
      SHARED ${fsunlinsolklu_SOURCES})
    set_target_properties(sundials_fsunlinsolklu_shared
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsolklu
      CLEAN_DIRECT_OUTPUT 1
      VERSION ${sunlinsollib_VERSION}
      SOVERSION ${sunlinsollib_SOVERSION})
    # depends on fnvecserial, fsunmatrixsparse, sunlinsolklu
    target_link_libraries(sundials_fsunlinsolklu_shared
      PUBLIC
      sundials_fnvecserial_shared
      sundials_fsunmatrixsparse_shared
      sundials_sunlinsolklu_shared)
    install(TARGETS sundials_fsunlinsolklu_shared
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

  #
  message(STATUS "Added SUNLINSOL_KLU F77 interface")

endif(F77_INTERFACE_ENABLE AND F77_FOUND)

# Add F90 module if F2003 interface is enabled
if(F2003_FOUND AND F2003_INTERFACE_ENABLE)
  add_subdirectory(fmod)
endif(F2003_FOUND AND F2003_INTERFACE_ENABLE)
