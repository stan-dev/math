
// Copyright 2022 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#if !defined(BOOST_UNORDERED_TEST_HELPERS_UNORDERED_HEADER)
#define BOOST_UNORDERED_TEST_HELPERS_UNORDERED_HEADER

// clang-format off
#include "prefix.hpp"
#ifdef BOOST_UNORDERED_FOA_TESTS
#include <boost/unordered/unordered_flat_set.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_node_map.hpp>
#include <boost/unordered/unordered_node_set.hpp>
#include <boost/unordered/detail/implementation.hpp>
#else
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#endif
#include "postfix.hpp"
// clang-format on

#if defined(BOOST_LIBSTDCXX_VERSION)
#if BOOST_LIBSTDCXX_VERSION < 60000
#define BOOST_UNORDERED_NO_INIT_TYPE_TESTS
#endif
#endif

#endif
