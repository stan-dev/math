# Copyright Hans Dembinski 2019
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

boost_test(SOURCES check_cmake_version.cpp ARGUMENTS ${PROJECT_VERSION}
  LIBRARIES Boost::core Boost::config)

# checks that b2 and cmake are in sync
add_test(NAME check_build_system COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/check_build_system.py)

# checks that all headers are included in odr test
add_test(NAME check_odr_test COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/check_odr_test.py)

# keep in sync with Jamfile, this should be automatized...
boost_test(TYPE compile-fail SOURCES axis_category_fail0.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES axis_category_fail1.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES axis_category_fail2.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES axis_integer_fail0.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES axis_integer_fail1.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES axis_integer_fail2.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES axis_integer_fail3.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES axis_integer_fail4.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES axis_regular_fail0.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES axis_regular_fail1.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES axis_variable_fail0.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES axis_variable_fail1.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES make_histogram_fail0.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES make_histogram_fail1.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES histogram_fail0.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES histogram_fail1.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES histogram_fail2.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES histogram_fail3.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE compile-fail SOURCES histogram_fail4.cpp
  LIBRARIES Boost::histogram
)
boost_test(TYPE run SOURCES accumulators_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES algorithm_project_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES algorithm_reduce_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES algorithm_sum_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES algorithm_empty_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES axis_category_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES axis_integer_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES axis_option_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES axis_regular_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES axis_size.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES axis_traits_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES axis_variable_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES axis_variant_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_accumulator_traits_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_argument_traits_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_args_type_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_axes_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_convert_integer_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_detect_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_limits_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_make_default_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_misc_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_large_int_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_iterator_adaptor_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_operators_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_relaxed_equal_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_replace_type_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_safe_comparison_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_static_if_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES detail_tuple_slice_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES histogram_custom_axis_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES histogram_dynamic_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES histogram_fill_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES histogram_growing_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES histogram_mixed_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES histogram_operators_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES histogram_ostream_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES histogram_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES indexed_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES storage_adaptor_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES unlimited_storage_test.cpp
  LIBRARIES Boost::histogram Boost::core)
boost_test(TYPE run SOURCES utility_test.cpp
  LIBRARIES Boost::histogram Boost::core)

if (cxx_std_17 IN_LIST CMAKE_CXX_COMPILE_FEATURES)
  boost_test(TYPE run SOURCES deduction_guides_test.cpp
    LIBRARIES Boost::histogram Boost::core)
  target_compile_features(BoostHistogram-deduction_guides_test_cpp PRIVATE cxx_std_17)
endif()

if (Threads_FOUND)
  boost_test(TYPE run SOURCES histogram_threaded_test.cpp
    LIBRARIES Boost::histogram Boost::core Threads::Threads)
  boost_test(TYPE run SOURCES storage_adaptor_threaded_test.cpp
    LIBRARIES Boost::histogram Boost::core Threads::Threads)
endif()

## No cmake support yet
# boost_test(TYPE link SOURCES odr_main_test.cpp odr_test.cpp
#  LIBRARIES Boost::histogram Boost::core Boost::serialization)
# boost_test(TYPE run SOURCES boost_accumulators_support_test.cpp
#  LIBRARIES Boost::histogram Boost::core Boost::accumulators)
# boost_test(TYPE run SOURCES boost_range_support_test.cpp
#  LIBRARIES Boost::histogram Boost::core Boost::range)
# boost_test(TYPE run SOURCES boost_units_support_test.cpp
#  LIBRARIES Boost::histogram Boost::core Boost::units)
# boost_test(TYPE run SOURCES detail_array_wrapper_serialization_test.cpp LIBRARIES Boost::histogram Boost::core Boost::serialization)
# boost_test(TYPE run SOURCES unlimited_storage_serialization_test.cpp LIBRARIES Boost::histogram Boost::core Boost::serialization)
# boost_test(TYPE run SOURCES storage_adaptor_serialization_test.cpp LIBRARIES Boost::histogram Boost::core Boost::serialization)
# boost_test(TYPE run SOURCES histogram_serialization_test.cpp LIBRARIES Boost::histogram Boost::core Boost::serialization)
# boost_test(TYPE run SOURCES axis_variant_serialization_test.cpp
#   LIBRARIES Boost::histogram Boost::core Boost::serialization)
# boost_test(TYPE run SOURCES accumulators_serialization_test.cpp
#   LIBRARIES Boost::histogram Boost::core Boost::serialization)
