# Copyright 2018-2020 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

include(BoostTest OPTIONAL RESULT_VARIABLE HAVE_BOOST_TEST)

if(NOT HAVE_BOOST_TEST)
  return()
endif()

macro(typeof_test s1)

  string(MAKE_C_IDENTIFIER ${s1} n1)

  boost_test(TYPE compile SOURCES ${s1} ${ARGN} COMPILE_DEFINITIONS BOOST_TYPEOF_NATIVE NAME ${n1}_native)
  boost_test(TYPE compile SOURCES ${s1} ${ARGN} COMPILE_DEFINITIONS BOOST_TYPEOF_EMULATION NAME ${n1}_emulation)

endmacro()

set(BOOST_TEST_LINK_LIBRARIES Boost::typeof Boost::core Boost::type_traits Boost::static_assert)

set(tests

data_member.cpp
function.cpp
function_binding.cpp
function_ptr.cpp
function_ptr_from_tpl.cpp
function_ref.cpp
member_function.cpp
modifiers.cpp
msvc_typeof_in_lambda.cpp
nested_typedef.cpp
noncopyable.cpp
std.cpp
template_dependent.cpp
template_enum.cpp
template_int.cpp
template_multiword.cpp
template_tpl.cpp
template_type.cpp
type.cpp
)

foreach(test IN LISTS tests)

  typeof_test(${test})

endforeach()

boost_test(TYPE run SOURCES odr1.cpp odr2.cpp COMPILE_DEFINITIONS BOOST_TYPEOF_NATIVE NAME odr_native)
boost_test(TYPE run SOURCES odr1.cpp odr2.cpp COMPILE_DEFINITIONS BOOST_TYPEOF_EMULATED NAME odr_emulated)

boost_test(TYPE run SOURCES odr_no_uns1.cpp odr_no_uns2.cpp COMPILE_DEFINITIONS BOOST_TYPEOF_EMULATED NAME odr_no_uns)
