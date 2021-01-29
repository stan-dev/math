# Copyright 2018, 2019 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

include(BoostTest OPTIONAL RESULT_VARIABLE HAVE_BOOST_TEST)

if(NOT HAVE_BOOST_TEST)
  return()
endif()

enable_language(C)

set(tests_common arithmetic array comparison control debug facilities list logical punctuation selection seq slot stringize tuple variadic isempty)
set(tests_c ${tests_common})
set(tests_cpp ${tests_common} iteration repetition quick)

set(BOOST_TEST_LINK_LIBRARIES Boost::preprocessor)

include_directories(../../..) # for `include <libs/preprocessor/...>` to work

foreach(test IN LISTS tests_c)

  boost_test(TYPE compile SOURCES ${test}.c)
  boost_test(TYPE compile NAME ${test}_c_nvm SOURCES ${test}.c COMPILE_DEFINITIONS BOOST_PP_VARIADICS=0)

endforeach()

foreach(test IN LISTS tests_cpp)

  boost_test(TYPE compile SOURCES ${test}.cpp)
  boost_test(TYPE compile NAME ${test}_cpp_nvm SOURCES ${test}.cpp COMPILE_DEFINITIONS BOOST_PP_VARIADICS=0)

endforeach()
