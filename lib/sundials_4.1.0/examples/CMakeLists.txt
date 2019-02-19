# ---------------------------------------------------------------
# Programmer:  David Gardner and Slaven Peles @ LLNL
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
# examples level CMakeLists.txt for SUNDIALS (for cmake build system)
# ---------------------------------------------------------------

# The examples will be linked with the library corresponding to the build type.
# Whenever building shared libraries, use them to link the examples.
IF(BUILD_SHARED_LIBS)
  SET(LINK_LIBRARY_TYPE "shared")
ELSE()
  SET(LINK_LIBRARY_TYPE "static")
ENDIF()

# Set variables used in generating CMake and Makefiles for examples
IF(EXAMPLES_INSTALL)

  SET(SHELL "sh")
  SET(prefix "${CMAKE_INSTALL_PREFIX}")
  SET(exec_prefix "${CMAKE_INSTALL_PREFIX}")
  SET(includedir "${prefix}/include")
  SET(libdir "${exec_prefix}/${CMAKE_INSTALL_LIBDIR}")
  SET(CPP "${CMAKE_C_COMPILER}")
  SET(CPPFLAGS "${CMAKE_C_FLAGS_RELEASE}")
  SET(CC "${CMAKE_C_COMPILER}")
  SET(CFLAGS "${CMAKE_C_FLAGS_RELEASE}")
  SET(LDFLAGS "${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
  LIST2STRING(EXTRA_LINK_LIBS LIBS)

  IF(CXX_FOUND)
    SET(CXX "${CMAKE_CXX_COMPILER}")
    SET(CXX_LNKR "${CMAKE_CXX_COMPILER}")
    SET(CXXFLAGS "${CMAKE_CXX_FLAGS_RELEASE}")
    SET(CXX_LDFLAGS "${CMAKE_CXX_FLAGS_RELEASE}")
    LIST2STRING(EXTRA_LINK_LIBS CXX_LIBS)
  ENDIF(CXX_FOUND)

  IF(F77_FOUND)
    SET(F77 "${CMAKE_Fortran_COMPILER}")
    SET(F77_LNKR "${CMAKE_Fortran_COMPILER}")
    SET(FFLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
    SET(F77_LDFLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
    LIST2STRING(EXTRA_LINK_LIBS F77_LIBS)
  ENDIF(F77_FOUND)

  IF(F90_FOUND)
    SET(F90 "${CMAKE_Fortran_COMPILER}")
    SET(F90_LNKR "${CMAKE_Fortran_COMPILER}")
    SET(F90FLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
    SET(F90_LDFLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
    LIST2STRING(EXTRA_LINK_LIBS F90_LIBS)
  ENDIF(F90_FOUND)

  IF(SUPERLUMT_FOUND)
    LIST2STRING(SUPERLUMT_LIBRARIES SUPERLUMT_LIBS)
    SET(SUPERLUMT_LIBS "${SUPERLUMT_LINKER_FLAGS} ${SUPERLUMT_LIBS}")
  ENDIF(SUPERLUMT_FOUND)

  IF(KLU_FOUND)
    LIST2STRING(KLU_LIBRARIES KLU_LIBS)
    SET(KLU_LIBS "${KLU_LINKER_FLAGS} ${KLU_LIBS}")
  ENDIF(KLU_FOUND)

  IF(BLAS_FOUND)
    LIST2STRING(BLAS_LIBRARIES BLAS_LIBS)
  ENDIF(BLAS_FOUND)

  IF(LAPACK_FOUND)
    LIST2STRING(LAPACK_LIBRARIES LAPACK_LIBS)
  ENDIF(LAPACK_FOUND)

  IF(MPI_C_FOUND)
    IF(MPI_C_COMPILER)
      SET(MPICC "${MPI_C_COMPILER}")
      SET(MPI_INC_DIR ".")
      SET(MPI_LIB_DIR ".")
      SET(MPI_LIBS "")
      SET(MPI_FLAGS "")
    ELSE()
      SET(MPICC "${CMAKE_C_COMPILER}")
      SET(MPI_INC_DIR "${MPI_INCLUDE_PATH}")
      SET(MPI_LIB_DIR ".")
      LIST2STRING(MPI_LIBRARIES MPI_LIBS)
    ENDIF()
    SET(HYPRE_INC_DIR "${HYPRE_INCLUDE_DIR}")
    SET(HYPRE_LIB_DIR "${HYPRE_LIBRARY_DIR}")
    SET(HYPRE_LIBS "${HYPRE_LIBRARIES}")
  ENDIF()

  IF(MPI_CXX_FOUND)
    IF(MPI_CXX_COMPILER)
      SET(MPICXX "${MPI_CXX_COMPILER}")
    ELSE()
      SET(MPICXX "${CMAKE_CXX_COMPILER}")
      LIST2STRING(MPI_LIBRARIES MPI_LIBS)
    ENDIF()
  ENDIF()

  IF(MPI_Fortran_FOUND)
    IF(MPI_Fortran_COMPILER)
      SET(MPIF77 "${MPI_Fortran_COMPILER}")
      SET(MPIF77_LNKR "${MPI_Fortran_COMPILER}")

      SET(MPIF90 "${MPI_Fortran_COMPILER}")
      SET(MPIF90_LNKR "${MPI_Fortran_COMPILER}")
    ELSE()
      SET(MPIF77 "${CMAKE_Fortran_COMPILER}")
      SET(MPIF77_LNKR "${CMAKE_Fortran_COMPILER}")
      SET(MPI_INC_DIR "${MPI_INCLUDE_PATH}")
      SET(MPI_LIB_DIR ".")
      LIST2STRING(MPI_LIBRARIES MPI_LIBS)

      SET(MPIF90 "${CMAKE_Fortran_COMPILER}")
      SET(MPIF90_LNKR "${CMAKE_Fortran_COMPILER}")
      LIST2STRING(MPI_LIBRARIES MPI_LIBS)
    ENDIF()
  ENDIF()

ENDIF(EXAMPLES_INSTALL)


#----------------------------------------
# Add specific examples
#----------------------------------------

# Add ARKode examples
IF(BUILD_ARKODE)
  ADD_SUBDIRECTORY(arkode)
ENDIF(BUILD_ARKODE)

# Add CVODE examples
IF(BUILD_CVODE)
  ADD_SUBDIRECTORY(cvode)
ENDIF(BUILD_CVODE)

# Add CVODES Examples
IF(BUILD_CVODES)
  ADD_SUBDIRECTORY(cvodes)
ENDIF(BUILD_CVODES)

# Add IDA examples
IF(BUILD_IDA)
  ADD_SUBDIRECTORY(ida)
ENDIF(BUILD_IDA)

# Add IDAS examples
IF(BUILD_IDAS)
  ADD_SUBDIRECTORY(idas)
ENDIF(BUILD_IDAS)

# Add KINSOL examples
IF(BUILD_KINSOL)
  ADD_SUBDIRECTORY(kinsol)
ENDIF(BUILD_KINSOL)

# Add CPODES examples
IF(BUILD_CPODES)
  ADD_SUBDIRECTORY(cpodes)
ENDIF(BUILD_CPODES)

# Add the nvector examples
ADD_SUBDIRECTORY(nvector)

# Add the sunmatrix examples
ADD_SUBDIRECTORY(sunmatrix)

# Add the sunlinearsolver examples
ADD_SUBDIRECTORY(sunlinsol)

# Add the sunnonlinearsolver examples
ADD_SUBDIRECTORY(sunnonlinsol)
