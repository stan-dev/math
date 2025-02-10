//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/error_code.hpp>

#include <boost/mysql/detail/intermediate_handler.hpp>

#include <boost/asio/associated_allocator.hpp>
#include <boost/asio/associated_cancellation_slot.hpp>
#include <boost/asio/associated_executor.hpp>
#include <boost/asio/associated_immediate_executor.hpp>
#include <boost/asio/async_result.hpp>
#include <boost/asio/bind_allocator.hpp>
#include <boost/asio/bind_cancellation_slot.hpp>
#include <boost/asio/bind_executor.hpp>
#include <boost/asio/bind_immediate_executor.hpp>
#include <boost/asio/cancellation_signal.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/test/unit_test.hpp>

#include "test_common/io_context_fixture.hpp"
#include "test_common/printing.hpp"
#include "test_common/tracker_executor.hpp"
#include "test_unit/custom_allocator.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
namespace asio = boost::asio;

BOOST_AUTO_TEST_SUITE(test_intermediate_handler)

//
// intermediate_handler. Verify that associated characteristics are propagated
//

// Initiation object that checks the expected values and calls the handler
struct initiate_check
{
    template <class Handler>
    void operator()(
        Handler&& handler,
        int expected_ex_id,
        int expected_immediate_ex_id,
        asio::cancellation_slot expected_slot
    )
    {
        // Associated executor propagated
        auto ex = asio::get_associated_executor(handler);
        BOOST_TEST(get_executor_id(ex) == expected_ex_id);

        // Immediate executor propagated
        auto immediate_ex = asio::get_associated_immediate_executor(handler, ex);
        BOOST_TEST(get_executor_id(immediate_ex) == expected_immediate_ex_id);

        // Cancellation slot propagated
        auto slot = asio::get_associated_cancellation_slot(handler);
        BOOST_TEST((slot == expected_slot));

        // Allocator propagated
        using alloc_t = decltype(asio::get_associated_allocator(handler));
        static_assert(std::is_same<alloc_t, custom_allocator<void>>::value, "");

        // Just call the handler
        handler(client_errc::wrong_num_params, 42);
    }
};

// Handler fn required by intermediate_handler
struct handler_fn
{
    int value;

    template <class Handler>
    void operator()(Handler&& handler, error_code ec, int v)
    {
        BOOST_TEST(ec == client_errc::wrong_num_params);
        BOOST_TEST(v == 42);
        handler(ec, v + value);
    }
};

// Initiates using an intermediate handler
template <class CompletionToken>
void async_check(
    int expected_ex_id,
    int expected_immediate_ex_id,
    asio::cancellation_slot expected_slot,
    CompletionToken&& token
)
{
    auto actual_handler = detail::make_intermediate_handler(
        handler_fn{3},
        std::forward<CompletionToken>(token)
    );
    asio::async_initiate<decltype(actual_handler), void(error_code, int)>(
        initiate_check{},
        actual_handler,
        expected_ex_id,
        expected_immediate_ex_id,
        expected_slot
    );
}

BOOST_FIXTURE_TEST_CASE(immediate_handler_propagates_properties, io_context_fixture)
{
    // Setup
    auto ex_result = create_tracker_executor(ctx.get_executor());
    auto immediate_ex_result = create_tracker_executor(ctx.get_executor());
    asio::cancellation_signal sig;
    custom_allocator<void> alloc;
    bool called = false;
    auto final_handler = [&](error_code ec, int value) {
        BOOST_TEST(ec == client_errc::wrong_num_params);
        BOOST_TEST(value == 45);
        called = true;
    };

    // Invoke the operation
    async_check(
        ex_result.executor_id,
        immediate_ex_result.executor_id,
        sig.slot(),
        asio::bind_executor(
            ex_result.ex,
            asio::bind_immediate_executor(
                immediate_ex_result.ex,
                asio::bind_cancellation_slot(
                    sig.slot(),
                    asio::bind_allocator(alloc, std::move(final_handler))
                )
            )
        )
    );

    // Sanity check
    BOOST_TEST(called);
}

BOOST_AUTO_TEST_SUITE_END()
