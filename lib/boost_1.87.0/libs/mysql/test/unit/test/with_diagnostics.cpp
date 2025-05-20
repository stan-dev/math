//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/with_diagnostics.hpp>

#include <boost/asio/async_result.hpp>
#include <boost/asio/bind_executor.hpp>
#include <boost/asio/cancel_after.hpp>
#include <boost/asio/consign.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/test/unit_test.hpp>

#include <exception>
#include <memory>
#include <type_traits>

#include "test_common/create_diagnostics.hpp"
#include "test_common/io_context_fixture.hpp"
#include "test_common/poll_until.hpp"
#include "test_common/printing.hpp"
#include "test_common/tracker_executor.hpp"
#include "test_unit/create_err.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/test_any_connection.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
namespace asio = boost::asio;

BOOST_AUTO_TEST_SUITE(test_with_diagnostics)

// Validate that an exception_ptr contains the expected values
void check_exception(std::exception_ptr exc, error_code expected_ec, const diagnostics& expected_diag)
{
    BOOST_CHECK_EXCEPTION(
        std::rethrow_exception(exc),
        error_with_diagnostics,
        [&](const error_with_diagnostics& e) {
            BOOST_TEST(e.code() == expected_ec);
            BOOST_TEST(e.get_diagnostics() == expected_diag);
            return true;
        }
    );
}

BOOST_FIXTURE_TEST_CASE(success, io_context_fixture)
{
    // Setup
    auto conn = create_test_any_connection(ctx);
    get_stream(conn).add_bytes(create_ok_frame(1, ok_builder().build()));
    bool called = false;

    // Execute
    conn.async_reset_connection(with_diagnostics([&](std::exception_ptr exc) {
        called = true;
        BOOST_TEST((exc == nullptr));
    }));
    poll_until(ctx, &called);
}

BOOST_FIXTURE_TEST_CASE(error, io_context_fixture)
{
    // Setup
    auto conn = create_test_any_connection(ctx);
    get_stream(conn).add_bytes(err_builder()
                                   .code(common_server_errc::er_no_such_user)
                                   .message("Invalid user")
                                   .seqnum(1)
                                   .build_frame());
    bool called = false;

    // Execute
    conn.async_reset_connection(with_diagnostics([&](std::exception_ptr exc) {
        called = true;
        check_exception(exc, common_server_errc::er_no_such_user, create_server_diag("Invalid user"));
    }));
    poll_until(ctx, &called);
}

struct nulldiag_initiation
{
    template <class Handler>
    void operator()(Handler&& handler, diagnostics* diag)
    {
        // with_diagnostics should have allocated a diagnostics object for us
        BOOST_TEST_REQUIRE(diag != nullptr);

        // set diagnostics
        *diag = create_server_diag("Invalid user");

        // call the handler
        std::move(handler)(common_server_errc::er_no_such_user);
    }
};

template <BOOST_ASIO_COMPLETION_TOKEN_FOR(void(error_code)) CompletionToken>
void async_nulldiag(CompletionToken&& token)
{
    diagnostics* diag = nullptr;
    asio::async_initiate<CompletionToken, void(error_code)>(nulldiag_initiation{}, token, diag);
}

BOOST_AUTO_TEST_CASE(diagnostics_null)
{
    // Setup
    bool called = false;

    // Execute
    async_nulldiag(with_diagnostics([&](std::exception_ptr exc) {
        called = true;
        check_exception(exc, common_server_errc::er_no_such_user, create_server_diag("Invalid user"));
    }));

    // Sanity check
    BOOST_TEST(called);
}

BOOST_FIXTURE_TEST_CASE(associated_properties, io_context_fixture)
{
    // Uses intermediate_handler internally, so we just perform a sanity check here
    // Setup
    auto ex_result = create_tracker_executor(ctx.get_executor());
    auto conn = create_test_any_connection(ctx);
    get_stream(conn).add_bytes(err_builder()
                                   .code(common_server_errc::er_no_such_user)
                                   .message("Invalid user")
                                   .seqnum(1)
                                   .build_frame());
    bool called = false;
    auto check_fn = [&](std::exception_ptr exc) {
        called = true;
        const int expected_stack[] = {ex_result.executor_id};
        BOOST_TEST(executor_stack() == expected_stack, boost::test_tools::per_element());
        check_exception(exc, common_server_errc::er_no_such_user, create_server_diag("Invalid user"));
    };

    // Call the op
    conn.async_reset_connection(with_diagnostics(asio::bind_executor(ex_result.ex, check_fn)));
    poll_until(ctx, &called);
}

// We correctly forward initiation args
struct test_initiation
{
    template <class Handler, class T1, class T2, class T3>
    void operator()(Handler&& handler, T1&& arg1, T2&& arg2, T3&& arg3, diagnostics*)
    {
        // T1 should be a non-const lvalue
        static_assert(std::is_same<T1, std::shared_ptr<int>&>::value, "");
        BOOST_TEST(arg1 != nullptr);

        // T2 should be a const lvalue
        static_assert(std::is_same<T2, const std::shared_ptr<int>&>::value, "");
        BOOST_TEST(arg2 != nullptr);

        // T3 should be a rvalue
        static_assert(std::is_same<T3, std::shared_ptr<int>>::value, "");
        BOOST_TEST(arg3 != nullptr);
        auto arg3_move = std::move(arg3);
        boost::ignore_unused(arg3_move);

        // Just call the handler
        std::move(handler)(error_code());
    }
};

template <BOOST_ASIO_COMPLETION_TOKEN_FOR(void(error_code)) CompletionToken>
void async_test(
    std::shared_ptr<int>& arg1,
    const std::shared_ptr<int>& arg2,
    std::shared_ptr<int>&& arg3,
    diagnostics& diag,
    CompletionToken&& token
)
{
    asio::async_initiate<CompletionToken, void(error_code)>(
        test_initiation{},
        token,
        arg1,
        arg2,
        std::move(arg3),
        &diag
    );
}

BOOST_AUTO_TEST_CASE(initiation_args_forwarding)
{
    // Setup
    auto arg1 = std::make_shared<int>(42);
    auto arg2 = arg1;
    auto arg3 = arg1;
    bool called = false;
    auto handler = [&](std::exception_ptr) { called = true; };
    diagnostics diag;

    // Call the operation
    async_test(arg1, arg2, std::move(arg3), diag, with_diagnostics(handler));

    // lvalues not moved, rvalue moved
    BOOST_TEST(arg1 != nullptr);
    BOOST_TEST(arg2 != nullptr);
    BOOST_TEST(arg3 == nullptr);
}

// Works fine if the token is a lvalue
BOOST_FIXTURE_TEST_CASE(token_lvalue, io_context_fixture)
{
    // Setup
    auto conn = create_test_any_connection(ctx);
    get_stream(conn).add_bytes(create_ok_frame(1, ok_builder().build()));
    bool called = false;
    auto token = with_diagnostics([&](std::exception_ptr) { called = true; });

    // Call the op
    conn.async_reset_connection(token);
    poll_until(ctx, &called);
}

BOOST_FIXTURE_TEST_CASE(token_const_lvalue, io_context_fixture)
{
    // Setup
    auto conn = create_test_any_connection(ctx);
    get_stream(conn).add_bytes(create_ok_frame(1, ok_builder().build()));
    bool called = false;
    const auto token = with_diagnostics([&](std::exception_ptr) { called = true; });

    // Call the op
    conn.async_reset_connection(token);
    poll_until(ctx, &called);
}

// with_diagnostics' initiation has the same executor as the initiation that gets passed,
// and thus with_diagnostics(asio::cancel_after(...)) works
BOOST_FIXTURE_TEST_CASE(initiation_propagates_executor, io_context_fixture)
{
    // Setup
    auto conn = create_test_any_connection(ctx);
    get_stream(conn).add_bytes(create_ok_frame(1, ok_builder().build()));
    bool called = false;
    const auto cb = [&](std::exception_ptr) { called = true; };

    // Call the op
    conn.async_reset_connection(with_diagnostics(asio::cancel_after(std::chrono::seconds(1), cb)));
    poll_until(ctx, &called);
}

// Edge case: if a diagnostics* gets passed as an argument
// to consign(), we don't mess things up
BOOST_FIXTURE_TEST_CASE(several_diagnostics_args, io_context_fixture)
{
    // Setup
    auto conn = create_test_any_connection(ctx);
    get_stream(conn).add_bytes(err_builder()
                                   .code(common_server_errc::er_no_such_user)
                                   .message("Invalid user")
                                   .seqnum(1)
                                   .build_frame());
    diagnostics other_diag;
    bool called = false;

    // Execute
    conn.async_reset_connection(asio::consign(
        with_diagnostics([&](std::exception_ptr exc) {
            called = true;
            check_exception(exc, common_server_errc::er_no_such_user, create_server_diag("Invalid user"));
        }),
        &other_diag
    ));
    poll_until(ctx, &called);

    // Unmodified
    BOOST_TEST(other_diag == diagnostics());
}

// with_diagnostics decays correctly
struct test_fn
{
    void operator()(error_code) {}
};
static_assert(
    std::is_same<decltype(with_diagnostics(std::declval<test_fn>())), with_diagnostics_t<test_fn>>::value,
    ""
);
static_assert(
    std::is_same<decltype(with_diagnostics(std::declval<test_fn&>())), with_diagnostics_t<test_fn>>::value,
    ""
);
static_assert(
    std::is_same<decltype(with_diagnostics(std::declval<const test_fn&>())), with_diagnostics_t<test_fn>>::
        value,
    ""
);
static_assert(
    std::is_same<decltype(with_diagnostics(std::declval<test_fn&&>())), with_diagnostics_t<test_fn>>::value,
    ""
);

// Applying with_diagnostics to an unkown signature is a pass-through
struct no_ec_initiation
{
    template <class Handler, class T1, class T2, class T3>
    void operator()(Handler&& handler, T1&& arg1, T2&& arg2, T3&& arg3)
    {
        // T1 should be a non-const lvalue
        static_assert(std::is_same<T1, std::shared_ptr<int>&>::value, "");
        BOOST_TEST(arg1 != nullptr);

        // T2 should be a const lvalue
        static_assert(std::is_same<T2, const std::shared_ptr<int>&>::value, "");
        BOOST_TEST(arg2 != nullptr);

        // T3 should be a rvalue
        static_assert(std::is_same<T3, std::shared_ptr<int>>::value, "");
        BOOST_TEST(arg3 != nullptr);
        auto arg3_move = std::move(arg3);
        boost::ignore_unused(arg3_move);

        // Just call the handler
        std::move(handler)(42);
    }
};

template <BOOST_ASIO_COMPLETION_TOKEN_FOR(void(int)) CompletionToken>
void async_no_ec(
    std::shared_ptr<int>& arg1,
    const std::shared_ptr<int>& arg2,
    std::shared_ptr<int>&& arg3,
    CompletionToken&& token
)
{
    asio::async_initiate<CompletionToken, void(int)>(no_ec_initiation{}, token, arg1, arg2, std::move(arg3));
}

BOOST_AUTO_TEST_CASE(signature_no_ec)
{
    // Setup
    auto arg1 = std::make_shared<int>(42);
    auto arg2 = arg1;
    auto arg3 = arg1;
    bool called = false;
    auto handler = [&](int val) {
        BOOST_TEST(val == 42);
        called = true;
    };

    // Call the operation
    async_no_ec(arg1, arg2, std::move(arg3), with_diagnostics(handler));

    // lvalues not moved, rvalue moved
    BOOST_TEST(arg1 != nullptr);
    BOOST_TEST(arg2 != nullptr);
    BOOST_TEST(arg3 == nullptr);
}

BOOST_AUTO_TEST_SUITE_END()
