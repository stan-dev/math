//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

// Test that header file is self-contained.
#include <boost/json/error.hpp>

#include <memory>

#include "test_suite.hpp"

namespace boost {
namespace json {

class error_test
{
public:
    void check(error e)
    {
        auto const ec = make_error_code(e);
        BOOST_TEST(ec.category().name() != nullptr);
        BOOST_TEST(! ec.message().empty());
        BOOST_TEST(ec.category().default_error_condition(
            static_cast<int>(e)).category() == ec.category());
    }

    void check(condition c, error e)
    {
        {
            auto const ec = make_error_code(e);
            BOOST_TEST(ec.category().name() != nullptr);
            BOOST_TEST(! ec.message().empty());
            BOOST_TEST(ec == c);
        }
        {
            auto ec = make_error_condition(c);
            BOOST_TEST(ec.category().name() != nullptr);
            BOOST_TEST(! ec.message().empty());
            BOOST_TEST(ec == c);
        }
    }

    void
    run()
    {
        check(condition::parse_error, error::syntax);
        check(condition::parse_error, error::extra_data);
        check(condition::parse_error, error::incomplete);
        check(condition::parse_error, error::exponent_overflow);
        check(condition::parse_error, error::too_deep);
        check(condition::parse_error, error::illegal_leading_surrogate);
        check(condition::parse_error, error::illegal_trailing_surrogate);
        check(condition::parse_error, error::expected_hex_digit);
        check(condition::parse_error, error::expected_utf16_escape);
        check(condition::parse_error, error::object_too_large);
        check(condition::parse_error, error::array_too_large);
        check(condition::parse_error, error::key_too_large);
        check(condition::parse_error, error::string_too_large);
        check(condition::parse_error, error::number_too_large);
        check(condition::parse_error, error::input_error);

        check(condition::pointer_parse_error, error::missing_slash);
        check(condition::pointer_parse_error, error::invalid_escape);

        check(condition::pointer_use_error, error::token_not_number);
        check(condition::pointer_use_error, error::value_is_scalar);
        check(condition::pointer_use_error, error::not_found);
        check(condition::pointer_use_error, error::token_overflow);
        check(condition::pointer_use_error, error::past_the_end);

        check(condition::conversion_error, error::not_number);
        check(condition::conversion_error, error::not_exact);
        check(condition::conversion_error, error::not_null);
        check(condition::conversion_error, error::not_bool);
        check(condition::conversion_error, error::not_array);
        check(condition::conversion_error, error::not_object);
        check(condition::conversion_error, error::not_string);
        check(condition::conversion_error, error::not_int64);
        check(condition::conversion_error, error::not_uint64);
        check(condition::conversion_error, error::not_double);
        check(condition::conversion_error, error::not_integer);
        check(condition::conversion_error, error::size_mismatch);
        check(condition::conversion_error, error::exhausted_variants);
        check(condition::conversion_error, error::unknown_name);

        check(condition::generic_error, error::exception);
        check(condition::generic_error, error::out_of_range);
        check(error::test_failure);

        // check std interop
        std::error_code const ec = error::syntax;
        BOOST_TEST(ec == error::syntax);
        BOOST_TEST(ec == condition::parse_error);
    }
};

TEST_SUITE(error_test, "boost.json.error");

} // namespace json
} // namespace boost
