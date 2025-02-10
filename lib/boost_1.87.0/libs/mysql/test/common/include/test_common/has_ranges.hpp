//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_HAS_RANGES_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_HAS_RANGES_HPP

#include <boost/config.hpp>

// libstdc++11 and below claim to support ranges, but basic piping fails to compile
#ifdef __cpp_lib_ranges
#if !defined(BOOST_LIBSTDCXX_VERSION) || BOOST_LIBSTDCXX_VERSION >= 120000
#define BOOST_MYSQL_HAS_RANGES
#endif
#endif

#endif
