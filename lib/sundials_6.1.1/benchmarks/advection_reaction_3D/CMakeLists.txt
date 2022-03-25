# ---------------------------------------------------------------
# Programmer(s):  Cody J. Balos @ LLNL
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

if(BUILD_ARKODE AND BUILD_CVODE AND BUILD_IDA)

  if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
    set(OTHER_LIBS OpenMP::OpenMP_CXX)
  endif()

  # ----------------------------------------------------------------------------
  # MPI only
  # ----------------------------------------------------------------------------

  add_executable(advection_reaction_3D
    advection_reaction_3D.cpp
    arkode_driver.cpp
    cvode_driver.cpp
    ida_driver.cpp
    rhs3D.hpp
    ParallelGrid.hpp
    backends.hpp)

  # ensure the linker language is reset to CXX
  set_target_properties(advection_reaction_3D PROPERTIES LINKER_LANGUAGE CXX)

  target_include_directories(advection_reaction_3D
    PRIVATE
    ${PROJECT_SOURCE_DIR}/utilities
    ${MPI_CXX_INCLUDE_DIRS})

  target_link_libraries(advection_reaction_3D
    PRIVATE
    sundials_arkode
    sundials_cvode
    sundials_ida
    sundials_nvecmpiplusx
    sundials_nvecserial
    RAJA
    ${MPI_CXX_LIBRARIES}
    ${OTHER_LIBS})

  install(TARGETS advection_reaction_3D
    DESTINATION "${CMAKE_INSTALL_BINDIR}/benchmarks/advection_reaction_3D")

  install(FILES README.md
    DESTINATION "${CMAKE_INSTALL_BINDIR}/benchmarks/advection_reaction_3D")

  # ----------------------------------------------------------------------------
  # MPI + CUDA
  # ----------------------------------------------------------------------------

  if(BUILD_NVECTOR_CUDA)

    set_source_files_properties(advection_reaction_3D.cpp
      PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(arkode_driver.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(cvode_driver.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(ida_driver.cpp PROPERTIES LANGUAGE CUDA)

    add_executable(advection_reaction_3D_mpicuda
      advection_reaction_3D.cpp
      arkode_driver.cpp
      cvode_driver.cpp
      ida_driver.cpp
      rhs3D.hpp
      ParallelGrid.hpp
      backends.hpp)

    # ensure the linker language is reset to CXX
    set_target_properties(advection_reaction_3D_mpicuda
      PROPERTIES LINKER_LANGUAGE CXX)

    target_include_directories(advection_reaction_3D_mpicuda
      PRIVATE
      ${PROJECT_SOURCE_DIR}/utilities
      ${MPI_CXX_INCLUDE_DIRS})

    target_link_libraries(advection_reaction_3D_mpicuda
      PRIVATE
      sundials_arkode
      sundials_cvode
      sundials_ida
      sundials_nvecmpiplusx
      sundials_nveccuda
      RAJA
      ${MPI_CXX_LIBRARIES}
      ${OTHER_LIBS})

    target_compile_definitions(advection_reaction_3D_mpicuda PRIVATE USE_CUDA_NVEC)

    install(TARGETS advection_reaction_3D_mpicuda
      DESTINATION "${CMAKE_INSTALL_BINDIR}/benchmarks/advection_reaction_3D")

  endif()

  # ----------------------------------------------------------------------------
  # MPI + HIP
  # ----------------------------------------------------------------------------

  if(BUILD_NVECTOR_HIP)

    add_executable(advection_reaction_3D_mpihip
      advection_reaction_3D.cpp
      arkode_driver.cpp
      cvode_driver.cpp
      ida_driver.cpp
      rhs3D.hpp
      ParallelGrid.hpp
      backends.hpp)

    target_include_directories(advection_reaction_3D_mpihip
      PRIVATE
      ${PROJECT_SOURCE_DIR}/utilities
      ${MPI_CXX_INCLUDE_DIRS})

    target_link_libraries(advection_reaction_3D_mpihip
      PRIVATE
      sundials_arkode
      sundials_cvode
      sundials_ida
      sundials_nvecmpiplusx
      sundials_nvechip
      RAJA
      hip::device
      ${MPI_CXX_LIBRARIES}
      ${OTHER_LIBS})

    target_compile_definitions(advection_reaction_3D_mpihip PRIVATE USE_HIP_NVEC)

    install(TARGETS advection_reaction_3D_mpihip
      DESTINATION "${CMAKE_INSTALL_BINDIR}/benchmarks/advection_reaction_3D")

  endif()

endif()
