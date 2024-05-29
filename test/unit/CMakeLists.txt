
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
