//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/pool_params.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/impl/internal/connection_pool/internal_pool_params.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/strand.hpp>
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <stdexcept>

using namespace boost::mysql;
namespace asio = boost::asio;

BOOST_AUTO_TEST_SUITE(test_pool_params)

BOOST_AUTO_TEST_CASE(invalid_params)
{
    struct
    {
        string_view name;
        void (*params_fn)(pool_params&);
        string_view expected_msg;
    } test_cases[] = {
        // clang-format off
        {
            "max_size 0",
            [](pool_params& p) { p.max_size = 0; },
            "pool_params::max_size must be greater than zero"
        },
        {
            "initial_size > max_size",
            [](pool_params& p) { p.max_size = 100; p.initial_size = 101; },
            "pool_params::max_size must be greater than pool_params::initial_size"
        },
        {
            "connect_timeout < 0",
            [](pool_params& p) { p.connect_timeout = std::chrono::seconds(-1); },
            "pool_params::connect_timeout must not be negative"
        },
        {
            "connect_timeout < 0 min",
            [](pool_params& p) { p.connect_timeout = (std::chrono::steady_clock::duration::min)(); },
            "pool_params::connect_timeout must not be negative"
        },
        {
            "retry_interval == 0",
            [](pool_params& p) { p.retry_interval = std::chrono::seconds(0); },
            "pool_params::retry_interval must be greater than zero"
        },
        {
            "retry_interval < 0",
            [](pool_params& p) { p.retry_interval = std::chrono::seconds(-1); },
            "pool_params::retry_interval must be greater than zero"
        },
        {
            "retry_interval < 0 min",
            [](pool_params& p) { p.retry_interval = (std::chrono::steady_clock::duration::min)(); },
            "pool_params::retry_interval must be greater than zero"
        },
        {
            "ping_interval < 0",
            [](pool_params& p) { p.ping_interval = std::chrono::seconds(-1); },
            "pool_params::ping_interval must not be negative"
        },
        {
            "ping_interval < 0 min",
            [](pool_params& p) { p.ping_interval = (std::chrono::steady_clock::duration::min)(); },
            "pool_params::ping_interval must not be negative"
        },
        {
            "ping_timeout < 0",
            [](pool_params& p) { p.ping_timeout = std::chrono::seconds(-1); },
            "pool_params::ping_timeout must not be negative"
        },
        {
            "ping_timeout < 0 min",
            [](pool_params& p) { p.ping_timeout = (std::chrono::steady_clock::duration::min)(); },
            "pool_params::ping_timeout must not be negative"
        },
        // clang-format on
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            pool_params params;
            tc.params_fn(params);
            auto matcher = [&tc](const std::invalid_argument& err) { return err.what() == tc.expected_msg; };
            BOOST_CHECK_EXCEPTION(detail::check_validity(params), std::invalid_argument, matcher);
        }
    }
}

BOOST_AUTO_TEST_CASE(valid_params)
{
    struct
    {
        string_view name;
        void (*params_fn)(pool_params&);
    } test_cases[] = {
        // clang-format off
        {
            "initial_size == 0",
            [](pool_params& p) { p.initial_size = 0; },
        },
        {
            "initial_size == max_size",
            [](pool_params& p) { p.max_size = 100; p.initial_size = 100; },
        },
        {
            "connect_timeout == 0",
            [](pool_params& p) { p.connect_timeout = std::chrono::seconds(0); },
        },
        {
            "connect_timeout == max",
            [](pool_params& p) { p.connect_timeout = (std::chrono::steady_clock::duration::max)(); },
        },
        {
            "retry_interval == max",
            [](pool_params& p) { p.retry_interval = (std::chrono::steady_clock::duration::max)(); },
        },
        {
            "ping_interval == 0",
            [](pool_params& p) { p.ping_interval = std::chrono::seconds(0); },
        },
        {
            "ping_interval == max",
            [](pool_params& p) { p.ping_interval = (std::chrono::steady_clock::duration::max)(); },
        },
        {
            "ping_timeout == 0",
            [](pool_params& p) { p.ping_timeout = std::chrono::seconds(0); },
        },
        {
            "ping_timeout == max",
            [](pool_params& p) { p.ping_timeout = (std::chrono::steady_clock::duration::max)(); },
        },
        {
            "thread_safe == true",
            [](pool_params& p) { p.thread_safe = true; },
        }
        // clang-format on
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            pool_params params;
            tc.params_fn(params);
            BOOST_CHECK_NO_THROW(detail::check_validity(params));
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
