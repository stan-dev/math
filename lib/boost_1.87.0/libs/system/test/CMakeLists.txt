# Copyright 2018-2020 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

include(BoostTest OPTIONAL RESULT_VARIABLE HAVE_BOOST_TEST)

if(NOT HAVE_BOOST_TEST)
  return()
endif()

macro(system_run s1)

  string(MAKE_C_IDENTIFIER ${s1} n1)

  boost_test(SOURCES ${s1} ${ARGN})
  boost_test(SOURCES ${s1} ${ARGN} COMPILE_DEFINITIONS BOOST_NO_ANSI_APIS NAME ${n1}_no_ansi)
  boost_test(SOURCES ${s1} ${ARGN} COMPILE_DEFINITIONS BOOST_SYSTEM_USE_UTF8 NAME ${n1}_utf8)
  boost_test(SOURCES ${s1} ${ARGN} COMPILE_DEFINITIONS BOOST_SYSTEM_DISABLE_THREADS NAME ${n1}_nthr)

endmacro()

macro(lib name macro)

  add_library(${name} EXCLUDE_FROM_ALL ${ARGN})
  target_compile_definitions(${name} PRIVATE $<$<BOOL:BUILD_SHARED_LIBS>:${macro}=1>)
  target_link_libraries(${name} PRIVATE Boost::system Boost::config)

endmacro()

set(BOOST_TEST_LINK_LIBRARIES Boost::system Boost::core Boost::static_assert)

boost_test(SOURCES quick.cpp)

system_run(error_code_test.cpp)
system_run(error_code_user_test.cpp)
system_run(system_error_test.cpp)

lib(boost_system_throw_test THROW_DYN_LINK throw_test.cpp)
boost_test(SOURCES dynamic_link_test.cpp LINK_LIBRARIES boost_system_throw_test)

system_run(initialization_test.cpp)
system_run(header_only_test.cpp)

boost_test(SOURCES config_test.cpp)

system_run(std_interop_test.cpp)
system_run(std_mismatch_test.cpp)

lib(boost_system_single_instance_lib1 SINGLE_INSTANCE_DYN_LINK single_instance_1.cpp)
lib(boost_system_single_instance_lib2 SINGLE_INSTANCE_DYN_LINK single_instance_2.cpp)

system_run(single_instance_test.cpp single_instance_1.cpp single_instance_2.cpp)
boost_test(SOURCES single_instance_test.cpp NAME single_instance_test_lib LINK_LIBRARIES boost_system_single_instance_lib1 boost_system_single_instance_lib2)

system_run(before_main_test.cpp)
system_run(constexpr_test.cpp)
system_run(win32_hresult_test.cpp)

system_run(error_category_test.cpp)
system_run(generic_category_test.cpp)
system_run(system_category_test.cpp)
system_run(after_main_test.cpp)
system_run(failed_test.cpp)
system_run(failed_constexpr_test.cpp)

boost_test(SOURCES warnings_test.cpp COMPILE_OPTIONS -Wall -Werror)

lib(boost_system_std_single_instance_lib1 STD_SINGLE_INSTANCE_DYN_LINK std_single_instance_1.cpp)
lib(boost_system_std_single_instance_lib2 STD_SINGLE_INSTANCE_DYN_LINK std_single_instance_2.cpp)

system_run(std_single_instance_test.cpp std_single_instance_1.cpp std_single_instance_2.cpp)
boost_test(SOURCES std_single_instance_test.cpp NAME std_single_instance_test_lib LINK_LIBRARIES boost_system_std_single_instance_lib1 boost_system_std_single_instance_lib2)

boost_test(TYPE run SOURCES is_error_code_enum_test.cpp)
boost_test(TYPE run SOURCES is_error_condition_enum_test.cpp)
boost_test(TYPE run SOURCES errc_test.cpp)
boost_test(TYPE run SOURCES errc_test2.cpp)
boost_test(TYPE run SOURCES error_category_test2.cpp)
boost_test(TYPE run SOURCES error_condition_test.cpp)
boost_test(TYPE run SOURCES error_condition_test2.cpp)
boost_test(TYPE run SOURCES generic_category_test2.cpp)
boost_test(TYPE run SOURCES generic_category_test3.cpp)
boost_test(TYPE run SOURCES system_category_test2.cpp)
boost_test(TYPE run SOURCES system_category_test3.cpp)

boost_test(TYPE run SOURCES windows_error_test.cpp)
boost_test(TYPE run SOURCES linux_error_test.cpp)

boost_test(TYPE link SOURCES errc_test3.cpp)

boost_test(TYPE run SOURCES snprintf_test.cpp)

boost_test(TYPE run SOURCES std_interop_test2.cpp)
boost_test(TYPE run SOURCES std_interop_test3.cpp)
boost_test(TYPE run SOURCES std_interop_test4.cpp)
boost_test(TYPE run SOURCES std_interop_test5.cpp)
boost_test(TYPE run SOURCES std_interop_test6.cpp)
boost_test(TYPE run SOURCES std_interop_test7.cpp)
boost_test(TYPE run SOURCES std_interop_test8.cpp)
boost_test(TYPE run SOURCES std_interop_test9.cpp)

boost_test(TYPE run SOURCES ec_location_test.cpp)

boost_test(TYPE run SOURCES error_condition_test3.cpp)
boost_test(TYPE run SOURCES error_code_test2.cpp)
boost_test(TYPE run SOURCES system_error_test2.cpp)
boost_test(TYPE run SOURCES std_interop_test10.cpp)

boost_test(TYPE run SOURCES ec_location_test2.cpp)
boost_test(TYPE run SOURCES ec_what_test.cpp)
boost_test(TYPE run SOURCES system_error_test3.cpp)

boost_test(TYPE run SOURCES std_interop_test11.cpp)

boost_test(TYPE run SOURCES ec_wstream_test.cpp)

boost_test(TYPE run SOURCES std_interop_test12.cpp)

boost_test(TYPE run SOURCES errc_test4.cpp)

boost_test(TYPE run SOURCES std_interop_test13.cpp)
boost_test(TYPE run SOURCES std_interop_test14.cpp)

boost_test(TYPE run SOURCES ec_location_test3.cpp)
boost_test(TYPE run SOURCES ec_location_test4.cpp)

boost_test(TYPE compile SOURCES constexpr_test2.cpp)

boost_test(TYPE run SOURCES error_code_test3.cpp)
boost_test(TYPE run SOURCES std_interop_test15.cpp)

boost_test(TYPE run SOURCES win32_generic_test.cpp)

boost_test(TYPE run SOURCES ec_hash_value_test.cpp)

boost_test(TYPE run SOURCES std_interop_test16.cpp)

# result

boost_test(TYPE run SOURCES result_default_construct.cpp)
boost_test(TYPE run SOURCES result_value_construct.cpp)
boost_test(TYPE run SOURCES result_error_construct.cpp)
boost_test(TYPE run SOURCES result_copy_construct.cpp)
boost_test(TYPE run SOURCES result_move_construct.cpp)
boost_test(TYPE run SOURCES result_copy_assign.cpp)
boost_test(TYPE run SOURCES result_move_assign.cpp)
boost_test(TYPE run SOURCES result_value_access.cpp)
boost_test(TYPE run SOURCES result_error_access.cpp)
boost_test(TYPE run SOURCES result_swap.cpp)
boost_test(TYPE run SOURCES result_eq.cpp)
boost_test(TYPE run SOURCES result_range_for.cpp)
boost_test(TYPE run SOURCES result_value_construct2.cpp)
boost_test(TYPE run SOURCES result_error_construct2.cpp)
boost_test(TYPE run SOURCES result_errc_construct.cpp)
boost_test(TYPE run SOURCES result_convert_construct.cpp)
boost_test(TYPE run SOURCES result_typedefs.cpp)
boost_test(TYPE run SOURCES result_value_construct3.cpp)
boost_test(TYPE run SOURCES result_error_construct3.cpp)
boost_test(TYPE run SOURCES result_emplace.cpp)
boost_test(TYPE run SOURCES result_error_construct4.cpp)
boost_test(TYPE run SOURCES result_value_construct4.cpp)
boost_test(TYPE run SOURCES result_value_construct5.cpp)
boost_test(TYPE run SOURCES result_error_move.cpp)
boost_test(TYPE run SOURCES result_value_construct6.cpp)
boost_test(TYPE run SOURCES result_value_construct7.cpp)
boost_test(TYPE run SOURCES result_error_construct5.cpp)
boost_test(TYPE run SOURCES result_or_value.cpp)
boost_test(TYPE compile-fail SOURCES result_or_value_fail.cpp)
boost_test(TYPE compile-fail SOURCES result_or_value_fail2.cpp)
boost_test(TYPE run SOURCES result_or_fn0v.cpp)
boost_test(TYPE run SOURCES result_or_fn0r.cpp)
boost_test(TYPE run SOURCES result_and_fn1v.cpp)
boost_test(TYPE run SOURCES result_and_fn1r.cpp)
boost_test(TYPE run SOURCES result_and_eq_fn1v.cpp)
boost_test(TYPE run SOURCES result_and_eq_fn1r.cpp)
boost_test(TYPE run SOURCES result_in_place_use.cpp)
boost_test(TYPE run SOURCES result_or_eq_value.cpp)
boost_test(TYPE run SOURCES result_or_eq_fn0v.cpp)
boost_test(TYPE run SOURCES result_or_eq_fn0r.cpp)
