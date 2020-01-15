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
# CMakeLists.txt file for the LapackBand SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_LAPACKBAND\n\")")

# Source files for the library
set(sunlinsollapackband_SOURCES sunlinsol_lapackband.c)

# Common SUNDIALS sources included in the library
set(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c)

# Exported header files
set(sunlinsollapackband_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_lapackband.h)

# Rules for building and installing the static library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Install the library
if(BUILD_STATIC_LIBS)

  add_library(sundials_sunlinsollapackband_static
    STATIC ${sunlinsollapackband_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsollapackband_static
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsollapackband
    CLEAN_DIRECT_OUTPUT 1)

  # depends on sunmatrixband and LAPACK
  target_link_libraries(sundials_sunlinsollapackband_static
    PUBLIC sundials_sunmatrixband_static ${LAPACK_LIBRARIES})

  target_compile_definitions(sundials_sunlinsollapackband_static
    PUBLIC -DBUILD_SUNDIALS_LIBRARY)

  install(TARGETS sundials_sunlinsollapackband_static
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the library
if(BUILD_SHARED_LIBS)

  add_library(sundials_sunlinsollapackband_shared
    SHARED ${sunlinsollapackband_SOURCES} ${shared_SOURCES})

  set_target_properties(sundials_sunlinsollapackband_shared
    PROPERTIES
    OUTPUT_NAME sundials_sunlinsollapackband
    CLEAN_DIRECT_OUTPUT 1
    VERSION ${sunlinsollib_VERSION}
    SOVERSION ${sunlinsollib_SOVERSION})

  # depends on sunmatrixband and LAPACK
  target_link_libraries(sundials_sunlinsollapackband_shared
    PUBLIC sundials_sunmatrixband_shared ${LAPACK_LIBRARIES})

  install(TARGETS sundials_sunlinsollapackband_shared
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_SHARED_LIBS)

# Install the header files
install(FILES ${sunlinsollapackband_HEADERS} DESTINATION include/sunlinsol)

#
message(STATUS "Added SUNLINSOL_LAPACKBAND module")

# If FCMIX is enabled, build and install the F77 library
if(F77_INTERFACE_ENABLE AND F77_FOUND)

  set(fsunlinsollapackband_SOURCES fsunlinsol_lapackband.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fsunlinsollapackband_static
      STATIC ${fsunlinsollapackband_SOURCES})
    set_target_properties(sundials_fsunlinsollapackband_static
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsollapackband
      CLEAN_DIRECT_OUTPUT 1)
    # depends on fnvecserial, fsunmatrixband, sunlinsollapackband
    target_link_libraries(sundials_fsunlinsollapackband_static
      PUBLIC
      sundials_fnvecserial_static
      sundials_fsunmatrixband_static
      sundials_sunlinsollapackband_static)
    install(TARGETS sundials_fsunlinsollapackband_static
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fsunlinsollapackband_shared
      SHARED ${fsunlinsollapackband_SOURCES})
    set_target_properties(sundials_fsunlinsollapackband_shared
      PROPERTIES
      OUTPUT_NAME sundials_fsunlinsollapackband
      CLEAN_DIRECT_OUTPUT 1
      VERSION ${sunlinsollib_VERSION}
      SOVERSION ${sunlinsollib_SOVERSION})
    # depends on fnvecserial, fsunmatrixband, sunlinsollapackband
    target_link_libraries(sundials_fsunlinsollapackband_shared
      PUBLIC
      sundials_fnvecserial_shared
      sundials_fsunmatrixband_shared
      sundials_sunlinsollapackband_shared)
    install(TARGETS sundials_fsunlinsollapackband_shared
      DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif(BUILD_SHARED_LIBS)

  #
  message(STATUS "Added SUNLINSOL_LAPACKBAND F77 interface")

endif(F77_INTERFACE_ENABLE AND F77_FOUND)
