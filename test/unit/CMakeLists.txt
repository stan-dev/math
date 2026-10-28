# Keep each autodiff layer's PCH small. The AD helper header includes mix.hpp
# and test_ad_matvar.hpp, so those need no separate entries in its PCH.
foreach(layer prim fwd rev mix)
  stan_add_pch(unit_${layer}_pch stan_test_options <gtest/gtest.h> <stan/math/${layer}.hpp>)
endforeach()
stan_add_pch(unit_ad_pch stan_test_options <test/unit/math/test_ad.hpp>)
add_gtest_grouped_test(PREFIX unit_math REPO_PATH test/unit/math
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/math" PCH unit_mix_pch
  LAYER_PCH_PREFIX unit AD_PCH unit_ad_pch)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  # Nonzero glibc allocations expose uninitialized CVODES quadrature buffers.
  set_tests_properties(test_unit_math_rev_functor_cos_ode_typed_test PROPERTIES
    ENVIRONMENT "MALLOC_PERTURB_=170")
endif()
# The timeout covers every parameterized case in the executable. Some
# Laplace suites exceed CTest's default 1500 seconds even with -O3.
get_property(laplace_tests DIRECTORY PROPERTY TESTS)
list(FILTER laplace_tests INCLUDE REGEX "^test_unit_math_laplace_")
if(laplace_tests)
  set_tests_properties(${laplace_tests} PROPERTIES TIMEOUT 3600)
endif()

# Umbrella and multiple-translation-unit checks must compile without a PCH.
# The public runner also accepts test/unit, so group them with test/unit/math.
set(unit_tests ${unit_math_selected_tests})
foreach(source math_c11_test.cpp math_include_test.cpp multiple_translation_units_test.cpp)
  get_filename_component(name "${source}" NAME_WE)
  stan_add_gtest(test_unit_${name} "${source}" "test/unit/${source}" "" unit)
  stan_matches_name_filters("test/unit/${source}" selected)
  if(selected)
    list(APPEND unit_tests test_unit_${name})
  endif()
endforeach()
add_library(stan_multiple_translation_units SHARED
  multiple_translation_units1.cpp multiple_translation_units2.cpp)
target_link_libraries(stan_multiple_translation_units PRIVATE stanmath)
target_link_libraries(test_unit_multiple_translation_units_test PRIVATE stan_multiple_translation_units)
if(unit_tests)
  stan_add_test_group(unit ${unit_tests})
  target_link_libraries(unit_all_tests PRIVATE stan_multiple_translation_units)
endif()
