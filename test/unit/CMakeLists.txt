
# Compile one test manually so that we can reuse precompile headers
add_executable(unit_pch
  ${CMAKE_CURRENT_SOURCE_DIR}/unit_pch.cpp)
target_include_directories(unit_pch PRIVATE ${CMAKE_SOURCE_DIR} ${boost_SOURCE_DIR})
target_precompile_headers(unit_pch PRIVATE
    "$<$<COMPILE_LANGUAGE:CXX>:<gtest/gtest.h$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:C>:<stddef.h$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/fwd.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/mix.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/rev.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:<stan/math/prim.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:test/unit/math/test_ad.hpp$<ANGLE-R>>"
    "$<$<COMPILE_LANGUAGE:CXX>:test/unit/math/test_ad_matvar.hpp$<ANGLE-R>>"
    )

target_link_libraries(unit_pch
  gtest_main benchmark::benchmark TBB::tbb
  Eigen3::Eigen sundials_kinsol sundials_cvodes sundials_ida
  sundials_idas sundials_nvecserial)



add_gtest_grouped_test(${CMAKE_CURRENT_SOURCE_DIR}/math unit_math unit_pch)

add_custom_target(unit_math_mix_subtests)
add_dependencies(unit_math_mix_subtests unit_math_mix_core unit_math_mix_err unit_math_mix_fun unit_math_mix_functor unit_math_mix_meta unit_math_mix_prob)
message("Adding unit_math_mix_subtests")

add_custom_target(unit_math_fwd_subtests)
add_dependencies(unit_math_fwd_subtests unit_math_fwd_core unit_math_fwd_err unit_math_fwd_fun unit_math_fwd_functor unit_math_fwd_meta unit_math_fwd_prob)
message("Adding unit_math_fwd_subtests")

add_custom_target(unit_math_prim_subtests)
add_dependencies(unit_math_prim_subtests unit_math_prim_core unit_math_prim_err unit_math_prim_fun unit_math_prim_functor unit_math_prim_meta unit_math_prim_prob)
message("Adding unit_math_prim_subtests")

add_custom_target(unit_math_rev_subtests)
add_dependencies(unit_math_rev_subtests unit_math_rev_core unit_math_rev_err unit_math_rev_fun unit_math_rev_functor unit_math_rev_meta unit_math_rev_prob)
message("Adding unit_math_rev_subtests")

