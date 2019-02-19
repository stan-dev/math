# CMakeLists.txt file for the FCVODE library

# Add variable fcvode_SOURCES with the sources for the FCVODE library
SET(fcvode_SOURCES
  fcvband.c
  fcvbbd.c
  fcvbp.c
  fcvdense.c
  fcvewt.c
  fcvjtimes.c
  fcvnulllinsol.c
  fcvnullmatrix.c
  fcvnullnonlinsol.c
  fcvode.c
  fcvpreco.c
  fcvroot.c
  fcvsparse.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the CVODE library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/nvector/serial/fnvector_serial.c
  )

# Add variable sunmatrix_SOURCES with the common SUNMatrix sources which will
# also be included in the CVODE library
SET(sunmatrix_SOURCES
  ${sundials_SOURCE_DIR}/src/sunmatrix/band/fsunmatrix_band.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/dense/fsunmatrix_dense.c
  ${sundials_SOURCE_DIR}/src/sunmatrix/sparse/fsunmatrix_sparse.c
  )

# Add variable sunlinsol_SOURCES with the common SUNLinearSolver sources which will
# also be included in the CVODE library
SET(sunlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunlinsol/band/fsunlinsol_band.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/dense/fsunlinsol_dense.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spbcgs/fsunlinsol_spbcgs.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spfgmr/fsunlinsol_spfgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/spgmr/fsunlinsol_spgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/sptfqmr/fsunlinsol_sptfqmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol/pcg/fsunlinsol_pcg.c
  )

IF(KLU_FOUND)
  LIST(APPEND sunlinsol_SOURCES
    ${sundials_SOURCE_DIR}/src/sunlinsol/klu/fsunlinsol_klu.c
    )
ENDIF()

IF(SUPERLUMT_FOUND)
  LIST(APPEND sunlinsol_SOURCES
    ${sundials_SOURCE_DIR}/src/sunlinsol/superlumt/fsunlinsol_superlumt.c
    )
ENDIF()

IF(LAPACK_FOUND)
  LIST(APPEND sunlinsol_SOURCES
    ${sundials_SOURCE_DIR}/src/sunlinsol/lapackband/fsunlinsol_lapackband.c
    ${sundials_SOURCE_DIR}/src/sunlinsol/lapackdense/fsunlinsol_lapackdense.c
    )
ENDIF()

# Add variable sunnonlinsol_SOURCES with the common SUNNonlinearSolver sources
# which will also be included in the CVODE library
SET(sunnonlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunnonlinsol/newton/fsunnonlinsol_newton.c
  ${sundials_SOURCE_DIR}/src/sunnonlinsol/fixedpoint/fsunnonlinsol_fixedpoint.c
  )

# Add source directories to include directories for access to
# implementation only header files (both for fcvode and cvode)
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(..)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Only build STATIC libraries (we cannot build shared libraries 
# for the FCMIX interfaces due to unresolved symbol errors 
# coming from inexistent user-provided functions)

# Add the build target for the FCVODE library
ADD_LIBRARY(sundials_fcvode_static STATIC
  ${fcvode_SOURCES} ${shared_SOURCES} ${sunmatrix_SOURCES} ${sunlinsol_SOURCES}
  ${sunnonlinsol_SOURCES})

# Set the library name and make sure it is not deleted
SET_TARGET_PROPERTIES(sundials_fcvode_static
  PROPERTIES OUTPUT_NAME sundials_fcvode CLEAN_DIRECT_OUTPUT 1)

# Install the FCVODE library
INSTALL(TARGETS sundials_fcvode_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

#
MESSAGE(STATUS "Added CVODE FCMIX module")
