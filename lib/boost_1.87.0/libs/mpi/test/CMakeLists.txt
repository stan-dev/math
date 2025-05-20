# Copyright 2023 Alain Miniussi
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

enable_testing()
set(CTEST_TEST_TIMEOUT 600)

function(add_mpi_tests_with_src name src)
  add_executable(${name} ${src})
  target_link_libraries(${name} PUBLIC Boost::mpi)
  set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  foreach(nproc ${ARGN})
    add_test(${name}-${nproc} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${nproc} ./${name})
  endforeach()
endfunction()

function(add_mpi_tests name)
  add_mpi_tests_with_src(${name} ${name}.cpp ${ARGN})
endfunction()

add_mpi_tests(test_all_reduce     1 2 7)
add_mpi_tests(test_block_nonblock 2)
add_mpi_tests(test_version 1)
add_mpi_tests_with_src(test_random_gather ../example/random_gather.cpp 2)
add_mpi_tests_with_src(test_random_scatter ../example/random_scatter.cpp 2)
add_mpi_tests_with_src(test_cartesian_communicator ../example/cartesian_communicator.cpp 24)
add_mpi_tests(test_cartesian_topology_init 1)
add_mpi_tests(test_broadcast_stl 2 )
add_mpi_tests(test_all_gather 1 2 11  )
add_mpi_tests(test_all_to_all 1 2 11  )
add_mpi_tests(test_broadcast 2 17 )
add_mpi_tests(test_gather 1 2 11  )
add_mpi_tests(test_is_mpi_op 1 )
add_mpi_tests(test_mt_level 1 )
 add_mpi_tests(test_mt_init 1 4 )
# # # Note: Microsoft MPI fails nonblocking_test on 1 processor
add_mpi_tests(test_nonblocking 2 11 24 )
add_mpi_tests(test_reduce 1 )
add_mpi_tests(test_ring 2 3 4 7 8 13 17 )
add_mpi_tests(test_sendrecv 1 4 7 48 )
add_mpi_tests(test_wait_any 1 4 7 20 )
add_mpi_tests(test_wait_all_vector 2 )
add_mpi_tests(test_wait_all_on_null 1 2 )
add_mpi_tests(test_scan 1 )
add_mpi_tests(test_scatter 1 )
# # # Note: Microsoft MPI fails all skeleton-content tests
add_mpi_tests(test_skeleton_content 2 3 4 7 8 13 17 )
add_mpi_tests(test_graph_topology 2 7 13 )
add_mpi_tests(test_cartesian_topology 24 )
add_mpi_tests(test_pointer 2 )
add_mpi_tests(test_groups 1 )
# # # tests that require -std=c++11
add_mpi_tests(test_sendrecv_vector 2 )
# # # Intel MPI 2018 and older are axtected to fail:
add_mpi_tests(test_non_blocking_any_source 2 17 )

