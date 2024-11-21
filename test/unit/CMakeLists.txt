
# Compile one test manually so that we can reuse precompile headers
add_executable(unit_pch
  ${CMAKE_CURRENT_SOURCE_DIR}/unit_pch.cpp)
target_include_directories(unit_pch PRIVATE ${CMAKE_SOURCE_DIR})
target_precompile_headers(unit_pch PRIVATE
"$<$<COMPILE_LANGUAGE:CXX>:<test/unit/math/test_ad.hpp$<ANGLE-R>>"
"$<$<COMPILE_LANGUAGE:CXX>:<test/unit/math/test_ad_matvar.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/mix.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/fwd.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/rev.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/prim.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<gtest/gtest.h$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:C>:<stddef.h$<ANGLE-R>>"
    )

target_link_libraries(unit_pch
  gtest_main benchmark::benchmark TBB::tbb
  Eigen3::Eigen sundials_kinsol sundials_cvodes sundials_ida
  sundials_idas sundials_nvecserial Boost::math Boost::numeric_odeint Boost::lexical_cast Boost::optional)

if (STAN_OPENCL)
target_include_directories(unit_pch PRIVATE ${CLCPP_SOURCE_DIR}/include)
target_link_libraries(unit_pch OpenCL::OpenCL)
endif()
if (STAN_MPI)
  target_link_libraries(unit_pch Boost::mpi MPI::MPI_CXX)
endif()

target_compile_options(unit_pch PRIVATE -O1)
add_gtest_grouped_test(${CMAKE_CURRENT_SOURCE_DIR}/math unit_math unit_pch)

message("Adding unit_math_opencl_subtests")
if (FALSE)


endif()


function(set_test_labels target label)
  get_property(SUBTESTS TARGET test_unit_math_prim_functor_subtests PROPERTY MANUALLY_ADDED_DEPENDENCIES)
  #loop over the subtests and and set the labels
  foreach(subtest ${SUBTESTS})
    set_property(TEST ${subtest} PROPERTY LABELS label)
  endforeach()
endfunction()

function(newset_test_labels target label)
  get_property(SUBTESTS TARGET ${target} PROPERTY MANUALLY_ADDED_DEPENDENCIES)
  # Loop over the subtests
  message(STATUS "{target} dependencies: ${SUBTESTS}")
  foreach(subtest ${SUBTESTS})
    string(FIND "${subtest}" "_test" INDEX REVERSE)
    message("Checking: ${subtest}")
    if(INDEX GREATER -1)
      message("${subtest} ends with the suffix.")
      # Set the labels for tests ending with "_test"
      set_property(TEST ${subtest} PROPERTY LABELS ${label})
    else()
      message(STATUS "RECURSE on ${subtest}")
      # Recursively call set_test_labels for other targets
      newset_test_labels(${subtest} ${label})
      message(STATUS "FINISH on ${subtest}")
    endif()
  endforeach()
endfunction()

add_custom_target(test_unit_math_opencl_tests)
add_dependencies(test_unit_math_opencl_tests
  test_unit_math_opencl_subtests test_unit_math_opencl_device_functions_subtests
  test_unit_math_opencl_kernel_generator_subtests
  test_unit_math_opencl_prim_subtests test_unit_math_opencl_rev_subtests)
newset_test_labels(test_unit_math_opencl_tests "unit_math_opencl")

message("Adding unit_math_mix_subtests")
add_custom_target(test_unit_math_mix_tests)
add_dependencies(test_unit_math_mix_tests
  test_unit_math_mix_core_subtests
  test_unit_math_mix_fun_subtests
  test_unit_math_mix_functor_subtests
  test_unit_math_mix_constraint_subtests
  test_unit_math_mix_meta_subtests
  test_unit_math_mix_prob_subtests)
newset_test_labels(test_unit_math_mix_tests "unit_math_mix")
newset_test_labels(test_unit_math_mix_fun_subtests "unit_math_mix_fun")

add_custom_target(test_unit_math_mix_nonfun_tests)
add_dependencies(test_unit_math_mix_nonfun_tests
  test_unit_math_mix_core_subtests
  test_unit_math_mix_functor_subtests
  test_unit_math_mix_meta_subtests
  test_unit_math_mix_prob_subtests)
newset_test_labels(test_unit_math_mix_nonfun_tests "unit_math_mix_nonfun")

message("Adding unit_math_fwd_nonfun_mix_subtests")
add_custom_target(test_unit_math_fwd_nonfun_mix_tests)
add_dependencies(test_unit_math_fwd_nonfun_mix_tests
  test_unit_math_fwd_core_subtests
  test_unit_math_fwd_fun_subtests
  test_unit_math_fwd_functor_subtests
  test_unit_math_fwd_meta_subtests
  test_unit_math_fwd_prob_subtests
  test_unit_math_mix_core_subtests
  test_unit_math_mix_functor_subtests
  test_unit_math_mix_meta_subtests
  test_unit_math_mix_prob_subtests)
newset_test_labels(test_unit_math_fwd_nonfun_mix_tests "unit_math_fwd_nonfun_mix")

message("Adding unit_math_fwd_subtests")
add_custom_target(test_unit_math_fwd_tests)
add_dependencies(test_unit_math_fwd_tests
  test_unit_math_fwd_core_subtests
  test_unit_math_fwd_fun_subtests
  test_unit_math_fwd_functor_subtests
  test_unit_math_fwd_meta_subtests
  test_unit_math_fwd_prob_subtests)

newset_test_labels(test_unit_math_fwd_tests "unit_math_fwd")

add_custom_target(test_unit_math_tests)
add_dependencies(test_unit_math_tests
  test_unit_math_subtests
  test_unit_math_fwd_subtests
  test_unit_math_fwd_core_subtests
  test_unit_math_fwd_fun_subtests
  test_unit_math_fwd_functor_subtests
  test_unit_math_fwd_meta_subtests
  test_unit_math_fwd_prob_subtests
  test_unit_math_memory_subtests
  test_unit_math_mix_subtests
  test_unit_math_mix_core_subtests
  test_unit_math_mix_fun_subtests
  test_unit_math_mix_functor_subtests
  test_unit_math_mix_meta_subtests
  test_unit_math_mix_prob_subtests
  test_unit_math_prim_core_subtests
  test_unit_math_prim_err_subtests
  test_unit_math_prim_fun_subtests
  test_unit_math_prim_functor_subtests
  test_unit_math_prim_meta_subtests
  test_unit_math_prim_prob_subtests
  test_unit_math_rev_subtests
  test_unit_math_rev_core_subtests
  test_unit_math_rev_err_subtests
  test_unit_math_rev_fun_subtests
  test_unit_math_rev_functor_subtests
  test_unit_math_rev_meta_subtests
  test_unit_math_rev_prob_subtests)
newset_test_labels(test_unit_math_tests "unit_math_subtest")



message("Adding unit_math_prim")
add_custom_target(test_unit_math_prim_tests)
add_dependencies(test_unit_math_prim_tests
  test_unit_math_prim_core_subtests
  test_unit_math_prim_err_subtests
  test_unit_math_prim_fun_subtests
  test_unit_math_prim_functor_subtests
  test_unit_math_prim_meta_subtests
  test_unit_math_prim_prob_subtests)
newset_test_labels(test_unit_math_prim_tests "unit_math_prim")

message("Adding unit_math_rev")
add_custom_target(test_unit_math_rev_tests)
add_dependencies(test_unit_math_rev_tests
  test_unit_math_rev_core_subtests
  test_unit_math_rev_err_subtests
  test_unit_math_rev_fun_subtests
  test_unit_math_rev_functor_subtests
  test_unit_math_rev_meta_subtests
  test_unit_math_rev_prob_subtests)
newset_test_labels(test_unit_math_rev_tests "unit_math_rev")


message("Adding test_unit_math_prim_rev_tests")
add_custom_target(test_unit_math_prim_rev_tests)
add_dependencies(test_unit_math_prim_rev_tests
  test_unit_math_prim_err_subtests
  test_unit_math_prim_fun_subtests
  test_unit_math_prim_functor_subtests
  test_unit_math_prim_meta_subtests
  test_unit_math_prim_prob_subtests
  test_unit_math_rev_core_subtests
  test_unit_math_rev_err_subtests
  test_unit_math_rev_fun_subtests
  test_unit_math_rev_functor_subtests
  test_unit_math_rev_meta_subtests
  test_unit_math_rev_prob_subtests)
newset_test_labels(test_unit_math_prim_rev_tests "unit_math_prim_rev")


message("Adding test_unit_math_mpi_subtests")
add_custom_target(test_unit_math_mpi_tests)
add_dependencies(test_unit_math_mpi_tests
  test_unit_math_prim_functor_subtests test_unit_math_rev_functor_subtests)
newset_test_labels(test_unit_math_mpi_tests "unit_math_mpi_subtest")


newset_test_labels(test_unit_math_prim_err_subtests "unit_math_prim_err")
