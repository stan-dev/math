# CMakeLists.txt file for the FCVODE library

# Add variable fcvode_SOURCES with the sources for the FCVODE library
SET(fcvode_SOURCES
  fcvband.c
  fcvbbd.c
  fcvbp.c
  fcvdense.c
  fcvewt.c
  fcvjtimes.c
  fcvode.c
  fcvpreco.c
  fcvroot.c
  )

IF(LAPACK_FOUND)
  SET(fcvode_BL_SOURCES fcvlapack.c fcvlapband.c fcvlapdense.c)
ELSE(LAPACK_FOUND)
  SET(fcvode_BL_SOURCES "")
ENDIF(LAPACK_FOUND)

IF(KLU_FOUND OR SUPERLUMT_FOUND)
  SET(fcvode_SPARSE_SOURCES fcvsparse.c)
ELSE(SPARSE_FOUND)
  SET(fcvode_SPARSE_SOURCES "")
ENDIF(KLU_FOUND OR SUPERLUMT_FOUND)

IF(KLU_FOUND)
  SET(fcvode_KLU_SOURCES fcvklu.c)
ELSE(KLU_FOUND)
  SET(fcvode_KLU_SOURCES "")
ENDIF(KLU_FOUND)

IF(SUPERLUMT_FOUND)
  SET(fcvode_SUPERLUMT_SOURCES fcvsuperlumt.c)
ELSE(SUPERLUMT_FOUND)
  SET(fcvode_SUPERLUMT_SOURCES "")
ENDIF(SUPERLUMT_FOUND)

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
ADD_LIBRARY(sundials_fcvode_static STATIC ${fcvode_SOURCES} ${fcvode_BL_SOURCES} ${fcvode_SPARSE_SOURCES} ${fcvode_KLU_SOURCES} ${fcvode_SUPERLUMT_SOURCES})

# Set the library name and make sure it is not deleted
SET_TARGET_PROPERTIES(sundials_fcvode_static
  PROPERTIES OUTPUT_NAME sundials_fcvode CLEAN_DIRECT_OUTPUT 1)

# Install the FCVODE library
INSTALL(TARGETS sundials_fcvode_static DESTINATION lib)

#
MESSAGE(STATUS "Added CVODE FCMIX module")
