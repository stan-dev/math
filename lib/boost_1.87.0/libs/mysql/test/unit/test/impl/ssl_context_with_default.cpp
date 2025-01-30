//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/impl/internal/ssl_context_with_default.hpp>

#include <boost/asio/ssl/context.hpp>
#include <boost/test/unit_test.hpp>

using namespace boost::mysql::detail;

BOOST_AUTO_TEST_SUITE(test_ssl_context_with_default)

BOOST_AUTO_TEST_CASE(no_external_context)
{
    // No external SSL context is passed
    ssl_context_with_default ctx_with_default(nullptr);

    // Nothing created on construction
    BOOST_TEST(ctx_with_default.get_ptr() == nullptr);

    // Calling get uses the default context singleton
    auto& ctx = ctx_with_default.get();
    auto handle = ctx.native_handle();
    BOOST_TEST(handle != nullptr);

    // Calling get again return the same context
    auto handle2 = ctx_with_default.get().native_handle();
    BOOST_TEST(handle == handle2);
}

BOOST_AUTO_TEST_CASE(external_context)
{
    // Create an external SSL context
    boost::asio::ssl::context ctx(boost::asio::ssl::context::tls_client);
    auto handle = ctx.native_handle();

    // Pass it to the object
    ssl_context_with_default ctx_with_default(&ctx);
    BOOST_TEST(ctx_with_default.get_ptr() == &ctx);

    // Calling get returns the passed context
    auto handle2 = ctx_with_default.get().native_handle();
    BOOST_TEST(handle == handle2);
}

BOOST_AUTO_TEST_SUITE_END()
