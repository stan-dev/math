# ---------------------------------------------------------------
# Programmer(s): David J. Gardner and Slaven Peles @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# examples level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# Set variables used in generating CMake and Makefiles for examples
if(EXAMPLES_INSTALL)

  set(SHELL "sh")
  set(prefix "${CMAKE_INSTALL_PREFIX}")
  set(exec_prefix "${CMAKE_INSTALL_PREFIX}")
  set(includedir "${prefix}/include")
  set(libdir "${exec_prefix}/${CMAKE_INSTALL_LIBDIR}")
  set(CPP "${_EXAMPLES_C_COMPILER}")
  set(CPPFLAGS "${CMAKE_C_FLAGS_RELEASE}")
  set(CC "${_EXAMPLES_C_COMPILER}")
  set(CFLAGS "${CMAKE_C_FLAGS_RELEASE}")
  set(LDFLAGS "${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
  list2string(EXTRA_LINK_LIBS LIBS)
  list2string(EXTRA_LINK_LIBS EXAMPLES_EXTRA_LIBS)

  if(CXX_FOUND)
    set(CXX "${_EXAMPLEX_CXX_COMPILER}")
    set(CXX_LNKR "${_EXAMPLEX_CXX_COMPILER}")
    set(CXXFLAGS "${CMAKE_CXX_FLAGS_RELEASE}")
    set(CXX_LDFLAGS "${CMAKE_CXX_FLAGS_RELEASE}")
    list2string(EXTRA_LINK_LIBS CXX_LIBS)
  endif()

  if(SUPERLUMT_FOUND)
    list2string(SUPERLUMT_LIBRARIES SUPERLUMT_LIBS)
    set(SUPERLUMT_LIBS "${SUPERLUMT_LINKER_FLAGS} ${SUPERLUMT_LIBS}")
  endif()

  if(KLU_FOUND)
    list2string(KLU_LIBRARIES KLU_LIBS)
    set(KLU_LIBS "${KLU_LINKER_FLAGS} ${KLU_LIBS}")
  endif()

  if(LAPACK_FOUND)
    list2string(LAPACK_LIBRARIES LAPACK_LIBS)
  endif()

  if(MPI_C_FOUND)
    if(MPI_C_COMPILER)
      set(MPICC "${MPI_C_COMPILER}")
      set(MPI_INC_DIR ".")
      set(MPI_LIB_DIR ".")
      set(MPI_LIBS "")
      set(MPI_FLAGS "")
    else()
      set(MPICC "${_EXAMPLES_C_COMPILER}")
      set(MPI_INC_DIR "${MPI_INCLUDE_PATH}")
      set(MPI_LIB_DIR ".")
      list2string(MPI_LIBRARIES MPI_LIBS)
    endif()
    set(HYPRE_INC_DIR "${HYPRE_INCLUDE_DIR}")
    set(HYPRE_LIB_DIR "${HYPRE_LIBRARY_DIR}")
    set(HYPRE_LIBS "${HYPRE_LIBRARIES}")
  endif()

  if(MPI_CXX_FOUND)
    if(MPI_CXX_COMPILER)
      set(MPICXX "${MPI_CXX_COMPILER}")
    else()
      set(MPICXX "${_EXAMPLEX_CXX_COMPILER}")
      list2string(MPI_LIBRARIES MPI_LIBS)
    endif()
  endif()

endif()


#----------------------------------------
# Add specific examples
#----------------------------------------

# Add ARKode examples
if(BUILD_ARKODE)
  add_subdirectory(arkode)
endif()

# Add CVODE examples
if(BUILD_CVODE)
  add_subdirectory(cvode)
endif()

# Add CVODES Examples
if(BUILD_CVODES)
  add_subdirectory(cvodes)
endif()

# Add IDA examples
if(BUILD_IDA)
  add_subdirectory(ida)
endif()

# Add IDAS examples
if(BUILD_IDAS)
  add_subdirectory(idas)
endif()

# Add KINSOL examples
if(BUILD_KINSOL)
  add_subdirectory(kinsol)
endif()

# Add CPODES examples
if(BUILD_CPODES)
  add_subdirectory(cpodes)
endif()

# Add the nvector examples
add_subdirectory(nvector)

# Add the sunmatrix examples
add_subdirectory(sunmatrix)

# Add the sunlinearsolver examples
add_subdirectory(sunlinsol)

# Add the sunnonlinearsolver examples
add_subdirectory(sunnonlinsol)
