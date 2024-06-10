
# Compile one test manually so that we can reuse precompile headers
add_executable(unit_pch
  ${CMAKE_CURRENT_SOURCE_DIR}/unit_pch.cpp)
target_include_directories(unit_pch PRIVATE ${CMAKE_SOURCE_DIR} ${boost_SOURCE_DIR})
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
  sundials_idas sundials_nvecserial)

if (STAN_OPENCL)
  target_link_libraries(unit_pch OpenCL::OpenCL)
  target_include_directories(unit_pch PRIVATE ${CLCPP_SOURCE_DIR}/include)
endif()
target_compile_options(unit_pch PRIVATE -O1)
add_gtest_grouped_test(${CMAKE_CURRENT_SOURCE_DIR}/math unit_math unit_pch)

message("Adding unit_math_opencl_subtests")
add_custom_target(unit_math_opencl_subtests)
add_dependencies(unit_math_opencl_subtests
test_unit_math_opencl test_unit_math_opencl_device_functions
test_unit_math_opencl_kernel_generator
test_unit_math_opencl_prim test_unit_math_opencl_rev)
set_property(TEST test_unit_math_opencl
  PROPERTY LABELS "unit_math_opencl")
set_property(TEST test_unit_math_opencl_device_functions
  PROPERTY LABELS "unit_math_opencl")
set_property(TEST test_unit_math_opencl_kernel_generator
  PROPERTY LABELS "unit_math_opencl")
set_property(TEST test_unit_math_opencl_prim
  PROPERTY LABELS "unit_math_opencl")
set_property(TEST test_unit_math_opencl_rev
  PROPERTY LABELS "unit_math_opencl")

message("Adding unit_math_mix_subtests")
add_custom_target(unit_math_mix_subtests)
add_dependencies(unit_math_mix_subtests test_unit_math_mix_core
  test_unit_math_mix_fun test_unit_math_mix_functor
  test_unit_math_mix_meta test_unit_math_mix_prob)
set_property(TEST test_unit_math_mix_core
  PROPERTY LABELS "unit_math_mix")
set_property(TEST test_unit_math_mix_fun
  PROPERTY LABELS "unit_math_mix")
set_property(TEST test_unit_math_mix_fun
  PROPERTY LABELS "unit_math_mix_fun")
set_property(TEST test_unit_math_mix_functor
  PROPERTY LABELS "unit_math_mix")
set_property(TEST test_unit_math_mix_meta
  PROPERTY LABELS "unit_math_mix")
set_property(TEST test_unit_math_mix_prob
  PROPERTY LABELS "unit_math_mix")

add_custom_target(unit_math_mix_nonfun_subtests)
add_dependencies(unit_math_mix_nonfun_subtests
  test_unit_math_mix_core test_unit_math_mix_functor
  test_unit_math_mix_meta test_unit_math_mix_prob)
set_property(TEST test_unit_math_mix_core
  PROPERTY LABELS "unit_math_mix_nonfun")
set_property(TEST test_unit_math_mix_functor
  PROPERTY LABELS "unit_math_mix_nonfun")
set_property(TEST test_unit_math_mix_meta
  PROPERTY LABELS "unit_math_mix_nonfun")
set_property(TEST test_unit_math_mix_prob
  PROPERTY LABELS "unit_math_mix_nonfun")

message("Adding unit_math_fwd_nonfun_mix_subtests")
add_custom_target(unit_math_fwd_nonfun_mix_subtests)
add_dependencies(unit_math_fwd_nonfun_mix_subtests
  test_unit_math_mix_core test_unit_math_mix_functor
  test_unit_math_mix_meta test_unit_math_mix_prob
  test_unit_math_fwd_core test_unit_math_fwd_fun
  test_unit_math_fwd_functor test_unit_math_fwd_meta
  test_unit_math_fwd_prob)
set_property(TEST test_unit_math_fwd_core
  PROPERTY LABELS "unit_math_fwd_nonfun_mix")
set_property(TEST test_unit_math_fwd_fun
  PROPERTY LABELS "unit_math_fwd_nonfun_mix")
set_property(TEST test_unit_math_fwd_functor
  PROPERTY LABELS "unit_math_fwd_nonfun_mix")
set_property(TEST test_unit_math_fwd_meta
  PROPERTY LABELS "unit_math_fwd_nonfun_mix")
set_property(TEST test_unit_math_fwd_prob
  PROPERTY LABELS "unit_math_fwd_nonfun_mix")
set_property(TEST test_unit_math_mix_core
  PROPERTY LABELS "unit_math_fwd_nonfun_mix")
set_property(TEST test_unit_math_mix_functor
  PROPERTY LABELS "unit_math_fwd_nonfun_mix")
set_property(TEST test_unit_math_mix_meta
  PROPERTY LABELS "unit_math_fwd_nonfun_mix")
set_property(TEST test_unit_math_mix_prob
  PROPERTY LABELS "unit_math_fwd_nonfun_mix")

message("Adding unit_math_fwd_subtests")
add_custom_target(unit_math_fwd_subtests)
add_dependencies(unit_math_fwd_subtests
  test_unit_math_fwd_core test_unit_math_fwd_fun
  test_unit_math_fwd_functor test_unit_math_fwd_meta
  test_unit_math_fwd_prob)
set_property(TEST test_unit_math_fwd_core
  PROPERTY LABELS "unit_math_fwd")
set_property(TEST test_unit_math_fwd_fun
  PROPERTY LABELS "unit_math_fwd")
set_property(TEST test_unit_math_fwd_functor
  PROPERTY LABELS "unit_math_fwd")
set_property(TEST test_unit_math_fwd_meta
  PROPERTY LABELS "unit_math_fwd")
set_property(TEST test_unit_math_fwd_prob
  PROPERTY LABELS "unit_math_fwd")

message("Adding unit_math_prim_subtests")
add_custom_target(unit_math_prim_subtests)
add_dependencies(unit_math_prim_subtests
  test_unit_math_prim_core test_unit_math_prim_err
  test_unit_math_prim_fun test_unit_math_prim_functor
  test_unit_math_prim_meta test_unit_math_prim_prob)
set_property(TEST test_unit_math_prim_core
  PROPERTY LABELS "unit_math_prim")
set_property(TEST test_unit_math_prim_err
  PROPERTY LABELS "unit_math_prim")
set_property(TEST test_unit_math_prim_fun
  PROPERTY LABELS "unit_math_prim")
set_property(TEST test_unit_math_prim_functor
  PROPERTY LABELS "unit_math_prim")
set_property(TEST test_unit_math_prim_meta
  PROPERTY LABELS "unit_math_prim")
set_property(TEST test_unit_math_prim_prob
  PROPERTY LABELS "unit_math_prim")

message("Adding unit_math_rev_subtests")
add_custom_target(unit_math_rev_subtests)
add_dependencies(unit_math_rev_subtests
  test_unit_math_rev_core test_unit_math_rev_err
  test_unit_math_rev_fun test_unit_math_rev_functor
  test_unit_math_rev_meta test_unit_math_rev_prob)
set_property(TEST test_unit_math_rev_core
  PROPERTY LABELS "unit_math_rev")
set_property(TEST test_unit_math_rev_err
  PROPERTY LABELS "unit_math_rev")
set_property(TEST test_unit_math_rev_fun
  PROPERTY LABELS "unit_math_rev")
set_property(TEST test_unit_math_rev_functor
  PROPERTY LABELS "unit_math_rev")
set_property(TEST test_unit_math_rev_meta
  PROPERTY LABELS "unit_math_rev")
set_property(TEST test_unit_math_rev_prob
  PROPERTY LABELS "unit_math_rev")

message("Adding unit_math_prim_rev_subtests")
add_custom_target(unit_math_prim_rev_subtests)
add_dependencies(unit_math_prim_rev_subtests
  test_unit_math_prim_core test_unit_math_prim_err
  test_unit_math_prim_fun test_unit_math_prim_functor
  test_unit_math_prim_meta test_unit_math_prim_prob
  test_unit_math_rev_core test_unit_math_rev_err
  test_unit_math_rev_fun test_unit_math_rev_functor
  test_unit_math_rev_meta test_unit_math_rev_prob)
set_property(TEST test_unit_math_prim_core
  PROPERTY LABELS "unit_math_prim_rev")
set_property(TEST test_unit_math_prim_err
  PROPERTY LABELS "unit_math_prim_rev")
set_property(TEST test_unit_math_prim_fun
  PROPERTY LABELS "unit_math_prim_rev")
set_property(TEST test_unit_math_prim_functor
  PROPERTY LABELS "unit_math_prim_rev")
set_property(TEST test_unit_math_prim_meta
  PROPERTY LABELS "unit_math_prim_rev")
set_property(TEST test_unit_math_prim_prob
  PROPERTY LABELS "unit_math_prim_rev")
set_property(TEST test_unit_math_rev_core
  PROPERTY LABELS "unit_math_prim_rev")
set_property(TEST test_unit_math_rev_err
  PROPERTY LABELS "unit_math_prim_rev")
set_property(TEST test_unit_math_rev_fun
  PROPERTY LABELS "unit_math_prim_rev")
set_property(TEST test_unit_math_rev_functor
  PROPERTY LABELS "unit_math_prim_rev")
set_property(TEST test_unit_math_rev_meta
  PROPERTY LABELS "unit_math_prim_rev")
set_property(TEST test_unit_math_rev_prob
  PROPERTY LABELS "unit_math_prim_rev")


  set_property(TEST test_unit_math PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_fwd PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_fwd_core PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_fwd_fun PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_fwd_functor PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_fwd_meta PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_fwd_prob PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_memory PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_mix PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_mix_core PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_mix_fun PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_mix_functor PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_mix_meta PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_mix_prob PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_prim_core PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_prim_err PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_prim_fun PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_prim_functor PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_prim_meta PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_prim_prob PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_rev PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_rev_core PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_rev_err PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_rev_fun PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_rev_functor PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_rev_meta PROPERTY LABELS "unit_math_subtest")
  set_property(TEST test_unit_math_rev_prob PROPERTY LABELS "unit_math_subtest")
