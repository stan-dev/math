# Copyright 2018-2020 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

include(BoostTest OPTIONAL RESULT_VARIABLE HAVE_BOOST_TEST)

if(NOT HAVE_BOOST_TEST)
  return()
endif()

boost_test(TYPE run SOURCES test.cpp test_impl.cpp LINK_LIBRARIES Boost::stacktrace Boost::core)
boost_test(TYPE run SOURCES test_noop.cpp test_impl.cpp LINK_LIBRARIES Boost::stacktrace_noop Boost::core)

boost_test(TYPE run SOURCES test_trivial.cpp LINK_LIBRARIES Boost::stacktrace Boost::core)
