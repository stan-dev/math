# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
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
# Defines the macro:
#
#     sundials_trycompile_execute(<EXECUTABLE> <CWD> <COMPILE_OK> <RUN_OK>
#                                 [COMPILE_OUTPUT variable]
#                                 [RUN_OUTPUT variable])
#
# This macro attempts to compile and then execute <CWD>/<EXECUTABLE>.
# The variable COMPILE_OK is TRUE if the source code compiles successfully.
# Otherwise COMPILE_OK is FALSE. The variable RUN_OK is TRUE if
# <CWD>/<EXECUTABLE> runs and returns zero. Otherwise it is FALSE.
# The optional COMPILE_OUTPUT variable is set to the generated output during
# compilation. It is useful for debugging compile failures. The option
# RUN_OUTPUT is set to the generated output during runtime.
# -----------------------------------------------------------------------------

macro(sundials_trycompile_execute EXECUTABLE CWD COMPILE_OK RUN_OK)

  set(options )
  set(oneValueArgs COMPILE_OUTPUT RUN_OUTPUT)
  set(multiValueArgs )
  set(COMPILE_OK FALSE)
  set(RUN_OK FALSE)
  set(COMPILE_OUTPUT )
  set(RUN_OUTPUT )

  cmake_parse_arguments(sundials_trycompile_execute "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  # compile the code and then try to run it
  try_compile(COMPILE_OK ${CWD} ${CWD} ${EXECUTABLE} OUTPUT_VARIABLE COMPILE_OUTPUT)
  if(COMPILE_OK)
    execute_process(COMMAND "./${EXECUTABLE}" WORKING_DIRECTORY ${CWD} RESULT_VARIABLE RUN_OK OUTPUT_VARIABLE RUN_OUTPUT)
    if(RUN_OK MATCHES "0")
      set(RUN_OK TRUE)
    endif()
  endif()

  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${CWD}/CMakeFiles)

  # set the optional outputs if used
  if(sundials_trycompile_execute_COMPILE_OUTPUT)
    set(${sundials_trycompile_execute_COMPILE_OUTPUT} ${COMPILE_OUTPUT})
  endif()
  if(sundials_trycompile_execute_RUN_OUTPUT)
    set(${sundials_trycompile_execute_RUN_OUTPUT} ${RUN_OUTPUT})
  endif()

  # unset variables we don't want to leak
  unset(COMPILE_OUTPUT)
  unset(RUN_OUTPUT)

endmacro()
