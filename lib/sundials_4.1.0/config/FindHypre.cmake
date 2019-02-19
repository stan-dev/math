# ---------------------------------------------------------------
# Programmer:  Slaven Peles @ LLNL, Jean Sexton @ SMU,
#              Eddy Banks @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# - Find hypre

#  HYPRE_INCLUDE_DIR = cached location of HYPRE.h
#  HYPRE_LIBRARY     = cached list of HYPRE library to link in

### Find include dir
find_path(temp_HYPRE_INCLUDE_DIR hypre.h ${HYPRE_INCLUDE_DIR})
if (temp_HYPRE_INCLUDE_DIR)
    set(HYPRE_INCLUDE_DIR ${temp_HYPRE_INCLUDE_DIR})
endif()
unset(temp_HYPRE_INCLUDE_DIR CACHE)
    
if (HYPRE_LIBRARY)
    # We have (or were given) HYPRE_LIBRARY - get path to use for any related libs
    get_filename_component(HYPRE_LIBRARY_DIR ${HYPRE_LIBRARY} PATH)

    # force CACHE update to show user DIR that will be used
    set(HYPRE_LIBRARY_DIR ${HYPRE_LIBRARY_DIR} CACHE PATH "" FORCE)
    
else ()
    # find library with user provided directory path
    set(HYPRE_LIBRARY_NAMES hypre HYPRE)
    find_library(HYPRE_LIBRARY 
      NAMES ${HYPRE_LIBRARY_NAMES}
      PATHS ${HYPRE_LIBRARY_DIR} NO_DEFAULT_PATH
      )
endif ()
mark_as_advanced(HYPRE_LIBRARY)

set(HYPRE_LIBRARIES ${HYPRE_LIBRARY})
