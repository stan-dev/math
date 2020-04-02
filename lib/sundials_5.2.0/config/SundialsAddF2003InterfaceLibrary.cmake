# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMake macro for adding F2003 interface libraries.
# Wraps the add_library command for adding Fortran 2003 interfaces.
#
#     sundials_add_f2003_interface_library(<target>
#         [STATIC | STATIC_OBJECT | SHARED | SHARED_OBJECT])
#
# It appends the library type to the Fortran_MODULE_DIRECTORY so
# that the .mod file generated during compilation of shared and
# static libraries do not collide. It also adds the .mod file
# directory to the target's includes. It also adds the C library
# as a library to link to.
#
# REQUIRES that <target> follows the SUNDIALS F2003 interface
# library naming convention i.e., if the C library/target is
# sundials_cvode_static then the <target> must be
# sundials_fcvode_mod_static.
#
# THE STATIC_OBJECT option generates a STATIC library target
# with the name <target> and an object library target with
# the name <target>_obj.
#
# THE SHARED_OBJECT option generates a SHARED library target
# with the name <target> and an object library target with
# the name <target>_obj.
# ---------------------------------------------------------------

macro(sundials_add_f2003_interface_library target)

  set(options STATIC STATIC_OBJECT SHARED SHARED_OBJECT)
  set(oneValueArgs )
  set(multiValueArgs )

  # parse keyword arguments/options
  cmake_parse_arguments(sundials_add_f2003_interface_library
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # source files are unparsed additional arguments
  set(sources ${sundials_add_f2003_interface_library_UNPARSED_ARGUMENTS})

  # determine the library type
  if(sundials_add_f2003_interface_library_STATIC OR sundials_add_f2003_interface_library_STATIC_OBJECT)
    set(libtype "STATIC")
  elseif(sundials_add_f2003_interface_library_SHARED OR sundials_add_f2003_interface_library_SHARED_OBJECT)
    set(libtype "SHARED")
  elseif(BUILD_SHARED_LIBS)
    set(libtype "SHARED")
  else()
    set(libtype "STATIC")
  endif()

  if(sundials_add_f2003_interface_library_STATIC_OBJECT OR sundials_add_f2003_interface_library_SHARED_OBJECT)
    # give the object library its own target name
    set(obj_target ${target}_obj)
    # create the target for the object library
    add_library(${obj_target} OBJECT ${sources})
    # now create the real library
    add_library(${target} ${libtype} $<TARGET_OBJECTS:${obj_target}>)
    # for shared libs, the objects must be PIC
    if("${libtype}" MATCHES "SHARED")
      set_target_properties(${obj_target} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    endif()
  else()
    # alias obj_target variable to target since there is no obj_target
    set(obj_target ${target})
    add_library(${target} ${libtype} ${sources})
  endif()

  # set target properties and target dependencies so that includes
  # and links get passed on when this target is used
  if(CMAKE_Fortran_MODULE_DIRECTORY)
    target_include_directories(${target}
      PUBLIC
        $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}_${libtype}>
        $<INSTALL_INTERFACE:${Fortran_INSTALL_MODDIR}>
      )
    set_target_properties(${obj_target} PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}_${libtype}")
    target_include_directories(${obj_target} PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY}_${libtype})
  else()
    target_include_directories(${target}
      PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${target}.dir>
        $<INSTALL_INTERFACE:${Fortran_INSTALL_MODDIR}>
      )
    set_target_properties(${obj_target} PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${target}.dir")
    target_include_directories(${obj_target} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${target}.dir)
  endif()

  # get the name of the C library which the fortran library interfaces to
  string(REPLACE "sundials_f" "sundials_" c_lib_name "${target}")
  string(REPLACE "_mod_" "_" c_lib_name "${c_lib_name}")

  # depend on the C library
  target_link_libraries(${target} PUBLIC ${c_lib_name})

endmacro(sundials_add_f2003_interface_library)
