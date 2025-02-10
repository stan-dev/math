//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_SOURCE_LOCATION_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_SOURCE_LOCATION_HPP

#include <boost/assert/source_location.hpp>

// boost::source_location triggers a bug on gcc < 8 when using PCHs
// BOOST_CURRENT_LOCATION complains about a redefinition of __PRETTY_FUNCTION__
// when used as default argument
#if defined(BOOST_GCC) && BOOST_GCC < 80000
#define BOOST_MYSQL_CURRENT_LOCATION ::boost::source_location(__FILE__, __LINE__, "")
#else
#define BOOST_MYSQL_CURRENT_LOCATION BOOST_CURRENT_LOCATION
#endif

#endif
