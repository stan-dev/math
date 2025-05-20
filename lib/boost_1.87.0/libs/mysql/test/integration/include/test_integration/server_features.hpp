//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SERVER_FEATURES_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SERVER_FEATURES_HPP

#include <boost/test/tree/decorator.hpp>

namespace boost {
namespace mysql {
namespace test {

// What does the CI server deployment support?
// We determine these from the command line arguments.
struct server_features
{
    // Is the server listening on a UNIX socket?
    bool unix_sockets{true};

    // Does the server support sha256 authentication methods?
    // Includes caching_sha2_password and sha256_password
    bool sha256{true};

    // Does the server support the dedicated JSON type?
    bool json_type{true};

    // Does the server support MySQL8+ specific regex error codes?
    bool regex_error_codes{true};

    // Does the server support MariaDB specific dup_query error codes?
    bool dup_query_error_codes{true};
};

using server_feature_t = bool server_features::*;

server_features get_server_features();

// Decorators to skip tests if the given features are not enabled.
// This is lazy because get_server_features requires access to getenv, not legal from global initializers.
unit_test::precondition run_if(server_feature_t feature);
unit_test::precondition run_if(server_feature_t feature1, server_feature_t feature2);

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
