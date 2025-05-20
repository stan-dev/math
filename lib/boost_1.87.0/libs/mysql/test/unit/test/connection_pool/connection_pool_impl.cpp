//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/character_set.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/connect_params.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/mysql_collations.hpp>
#include <boost/mysql/pipeline.hpp>
#include <boost/mysql/pool_params.hpp>
#include <boost/mysql/ssl_mode.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/access.hpp>
#include <boost/mysql/detail/pipeline.hpp>

#include <boost/mysql/impl/internal/connection_pool/connection_node.hpp>
#include <boost/mysql/impl/internal/connection_pool/connection_pool_impl.hpp>
#include <boost/mysql/impl/internal/connection_pool/internal_pool_params.hpp>
#include <boost/mysql/impl/internal/connection_pool/sansio_connection_node.hpp>

#include <boost/asio/any_completion_handler.hpp>
#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/async_result.hpp>
#include <boost/asio/bind_cancellation_slot.hpp>
#include <boost/asio/bind_executor.hpp>
#include <boost/asio/bind_immediate_executor.hpp>
#include <boost/asio/cancellation_signal.hpp>
#include <boost/asio/cancellation_type.hpp>
#include <boost/asio/compose.hpp>
#include <boost/asio/coroutine.hpp>
#include <boost/asio/dispatch.hpp>
#include <boost/asio/error.hpp>
#include <boost/asio/experimental/channel.hpp>
#include <boost/asio/experimental/channel_error.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/post.hpp>
#include <boost/asio/ssl/context.hpp>
#include <boost/asio/steady_timer.hpp>
#include <boost/asio/strand.hpp>
#include <boost/assert/source_location.hpp>
#include <boost/core/span.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <cstddef>
#include <memory>
#include <ostream>
#include <utility>

#include "test_common/create_diagnostics.hpp"
#include "test_common/io_context_fixture.hpp"
#include "test_common/network_result.hpp"
#include "test_common/poll_until.hpp"
#include "test_common/printing.hpp"
#include "test_common/source_location.hpp"
#include "test_common/tracker_executor.hpp"
#include "test_unit/mock_timer.hpp"
#include "test_unit/printing.hpp"

// These tests rely on channels, which are not compatible with this
// See https://github.com/chriskohlhoff/asio/issues/1398
#ifndef BOOST_ASIO_USE_TS_EXECUTOR_AS_DEFAULT

using namespace boost::mysql;
using namespace boost::mysql::test;
namespace asio = boost::asio;
using boost::test_tools::per_element;
using detail::connection_status;
using std::chrono::steady_clock;
namespace data = boost::unit_test::data;

/**
 * These tests verify step-by-step that all interactions between
 * elements in connection_pool work as intended. They use the templating
 * on IoTraits to mock out I/O objects. This allows for fast and reliable tests.
 * Boost.Asio lacks testing infrastructure, so we've coded some here.
 * These are complex tests.
 */

BOOST_AUTO_TEST_SUITE(test_connection_pool_impl)

enum class fn_type
{
    connect,
    pipeline,
    ping,
};

std::ostream& operator<<(std::ostream& os, fn_type t)
{
    switch (t)
    {
    case fn_type::connect: return os << "fn_type::connect";
    case fn_type::pipeline: return os << "fn_type::pipeline";
    case fn_type::ping: return os << "fn_type::ping";
    default: return os << "<unknown fn_type>";
    }
}

// A mock for mysql::any_connection. This allows us to control
// when and how operations like async_connect or async_ping complete,
// make assertions on the passed parameters, and force error conditions.
// By default, mocked operations stay outstanding until they're acknowledged
// by the test by calling mock_connection::step. step will wait until the appropriate
// mocked function is called and will make the outstanding operation complete with
// the passed error_code/diagnostics. All this synchronization uses channels.
class mock_connection
{
    struct impl_t
    {
        asio::experimental::channel<void(error_code, fn_type)> to_test_chan_;
        asio::experimental::channel<void(error_code, diagnostics)> from_test_chan_;

        // Transforms the channel-specific cancel code into asio::error::operation_aborted,
        // which is returned by connections when operations get cancelled
        static error_code transform_ec(error_code input)
        {
            return input == asio::experimental::channel_errc::channel_cancelled
                       ? asio::error::operation_aborted
                       : input;
        }

        // Code shared between all mocked ops
        struct mocked_op
        {
            impl_t& obj;
            fn_type op_type;
            diagnostics* diag;

            template <class Self>
            void operator()(Self& self)
            {
                // Notify the test that we're about to do op_type
                obj.to_test_chan_.async_send(error_code(), op_type, std::move(self));
            }

            template <class Self>
            void operator()(Self& self, error_code ec)
            {
                if (ec)
                {
                    // We were cancelled
                    self.complete(transform_ec(ec));
                }
                else
                {
                    // Read from the test what we should return
                    obj.from_test_chan_.async_receive(std::move(self));
                }
            }

            template <class Self>
            void operator()(Self& self, error_code ec, diagnostics recv_diag)
            {
                // Done
                if (diag)
                    *diag = std::move(recv_diag);
                self.complete(transform_ec(ec));
            }
        };

        template <class CompletionToken>
        auto op_impl(fn_type op_type, diagnostics* diag, CompletionToken&& token)
            -> decltype(asio::async_compose<CompletionToken, void(error_code)>(
                std::declval<mocked_op>(),
                token,
                asio::any_io_executor()
            ))
        {
            return asio::async_compose<CompletionToken, void(error_code)>(
                mocked_op{*this, op_type, diag},
                token,
                to_test_chan_.get_executor()
            );
        }

    } impl_;

    struct step_op
    {
        impl_t& obj;
        fn_type expected_op_type;
        error_code op_ec;
        diagnostics op_diag;

        template <class Self>
        void operator()(Self& self)
        {
            // Wait until the code under test performs the operation we want
            obj.to_test_chan_.async_receive(std::move(self));
        }

        template <class Self>
        void operator()(Self& self, error_code ec, fn_type actual_op_type)
        {
            // Verify it was actually what we wanted
            BOOST_TEST(ec == error_code());
            BOOST_TEST(actual_op_type == expected_op_type);

            // Tell the operation what its result should be
            obj.from_test_chan_.async_send(op_ec, std::move(op_diag), std::move(self));
        }

        template <class Self>
        void operator()(Self& self, error_code ec)
        {
            // Done
            BOOST_TEST(ec == error_code());
            self.complete(error_code());
        }
    };

    static void check_stages(const pipeline_request& req)
    {
        const detail::pipeline_request_stage expected_stages[] = {
            {detail::pipeline_stage_kind::reset_connection,  1, {}             },
            {detail::pipeline_stage_kind::set_character_set, 1, utf8mb4_charset},
        };
        BOOST_TEST(detail::access::get_impl(req).stages_ == expected_stages, per_element());
    }

    static void set_response(std::vector<stage_response>& res)
    {
        // Response should have two items, set to empty errors
        res.resize(2);
        for (auto& item : res)
            detail::access::get_impl(item).emplace_error();
    }

public:
    boost::mysql::any_connection_params ctor_params;
    boost::mysql::connect_params last_connect_params;

    mock_connection(asio::any_io_executor ex, boost::mysql::any_connection_params ctor_params)
        : impl_{ex, ex}, ctor_params(ctor_params)
    {
    }

    template <class CompletionToken>
    auto async_connect(const connect_params& params, diagnostics& diag, CompletionToken&& token)
        -> decltype(impl_.op_impl(fn_type::connect, &diag, std::forward<CompletionToken>(token)))
    {
        last_connect_params = params;
        return impl_.op_impl(fn_type::connect, &diag, std::forward<CompletionToken>(token));
    }

    template <class CompletionToken>
    auto async_ping(CompletionToken&& token
    ) -> decltype(impl_.op_impl(fn_type::ping, nullptr, std::forward<CompletionToken>(token)))
    {
        return impl_.op_impl(fn_type::ping, nullptr, std::forward<CompletionToken>(token));
    }

    template <class CompletionToken>
    auto async_run_pipeline(
        const pipeline_request& req,
        std::vector<stage_response>& res,
        CompletionToken&& token
    ) -> decltype(impl_.op_impl(fn_type::pipeline, nullptr, std::forward<CompletionToken>(token)))
    {
        check_stages(req);
        set_response(res);  // This should technically happen after initiation, but is enough for these tests
        return impl_.op_impl(fn_type::pipeline, nullptr, std::forward<CompletionToken>(token));
    }

    template <class CompletionToken>
    auto step(fn_type expected_op_type, CompletionToken&& token, error_code ec = {}, diagnostics diag = {})
        -> decltype(asio::async_compose<CompletionToken, void(error_code)>(
            step_op{impl_, expected_op_type, ec, std::move(diag)},
            token,
            impl_.from_test_chan_
        ))
    {
        return asio::async_compose<CompletionToken, void(error_code)>(
            step_op{impl_, expected_op_type, ec, std::move(diag)},
            token,
            impl_.from_test_chan_
        );
    }
};

struct mock_pooled_connection;
using mock_node = detail::basic_connection_node<mock_connection, mock_clock>;
using mock_pool = detail::basic_pool_impl<mock_connection, mock_clock, mock_pooled_connection>;

// Mock for pooled_connection
struct mock_pooled_connection
{
    std::shared_ptr<mock_pool> pool;
    mock_node* node{};

    mock_pooled_connection() = default;
    mock_pooled_connection(mock_node& node, std::shared_ptr<mock_pool> pool) noexcept
        : pool(std::move(pool)), node(&node)
    {
    }
};

// Helper class to launch an async_get_connection and then wait for it
// and validate its results
class get_connection_task
{
    struct impl_t
    {
        bool called{};
        mock_pool& pool;
        mock_node* actual_node{};
        mock_pool* actual_pool{};
        error_code actual_ec;
        bool was_immediate{};  // was the completion immediate?
        asio::cancellation_signal sig;

        impl_t(mock_pool& p) : pool(p) {}
    };

    struct handler
    {
        tracker_executor_result ex_result;
        tracker_executor_result immediate_ex_result;
        std::shared_ptr<impl_t> state;
        bool bind_slot;

        handler(std::shared_ptr<impl_t> st, bool bind_slot)
            : ex_result(create_tracker_executor(st->pool.get_executor())),
              immediate_ex_result(create_tracker_executor(st->pool.get_executor())),
              state(std::move(st)),
              bind_slot(bind_slot)
        {
        }

        using executor_type = asio::any_io_executor;
        executor_type get_executor() const { return ex_result.ex; }

        using immediate_executor_type = asio::any_io_executor;
        immediate_executor_type get_immediate_executor() const { return immediate_ex_result.ex; }

        using cancellation_slot_type = asio::cancellation_slot;
        cancellation_slot_type get_cancellation_slot() const
        {
            return bind_slot ? state->sig.slot() : asio::cancellation_slot();
        }

        void operator()(error_code ec, mock_pooled_connection conn)
        {
            // Check executor
            bool was_immediate = is_initiation_function();
            if (was_immediate)
            {
                // An immediate completion dispatches to the immediate executor,
                // then to the token's executor
                const int expected_stack[] = {immediate_ex_result.executor_id, ex_result.executor_id};
                BOOST_TEST(executor_stack() == expected_stack, per_element());
            }
            else
            {
                const int expected_stack[] = {ex_result.executor_id};
                BOOST_TEST(executor_stack() == expected_stack, per_element());
            }

            // If the pool was thread-safe, the callback should never be happening from
            // the pool's strand
            if (state->pool.params().thread_safe)
            {
                BOOST_TEST(!state->pool.strand().running_in_this_thread());
            }

            state->was_immediate = was_immediate;
            state->actual_node = conn.node;
            state->actual_pool = conn.pool.get();
            state->actual_ec = ec;
            state->called = true;
        }
    };

    std::shared_ptr<impl_t> impl_;

    void wait_impl(
        mock_node* expected_node,
        error_code expected_ec,
        bool expect_immediate,
        boost::source_location loc = BOOST_MYSQL_CURRENT_LOCATION
    )
    {
        auto& ctx = static_cast<asio::io_context&>(impl_->pool.get_executor().context());
        poll_until(ctx, &impl_->called, loc);
        BOOST_TEST_CONTEXT("Called from " << loc)
        {
            auto* expected_pool = expected_ec ? nullptr : &impl_->pool;
            BOOST_TEST(impl_->actual_ec == expected_ec);
            BOOST_TEST(impl_->actual_pool == expected_pool);
            BOOST_TEST(impl_->actual_node == expected_node);
            BOOST_TEST(impl_->was_immediate == expect_immediate);
        }
    }

public:
    get_connection_task(mock_pool& pool, diagnostics* diag, bool bind_slot = true)
        : impl_(std::make_shared<impl_t>(pool))
    {
        // Call the initiating function.
        // We need to dispatch it so the initiation runs in the io_context, and we can see immediate
        // completions
        auto impl = impl_;
        asio::dispatch(asio::bind_executor(pool.get_executor(), [impl, &pool, diag, bind_slot]() {
            // Mark that we're calling an initiating function
            initiation_guard guard;

            // Initiate
            pool.async_get_connection(diag, handler{impl, bind_slot});
        }));
    }

    void wait(
        mock_node& expected_node,
        bool expect_immediate,
        boost::source_location loc = BOOST_MYSQL_CURRENT_LOCATION
    )
    {
        wait_impl(&expected_node, error_code(), expect_immediate, loc);
    }

    void wait(
        error_code expected_ec,
        bool expect_immediate,
        boost::source_location loc = BOOST_MYSQL_CURRENT_LOCATION
    )
    {
        wait_impl(nullptr, expected_ec, expect_immediate, loc);
    }

    void cancel(asio::cancellation_type_t type = asio::cancellation_type_t::terminal)
    {
        impl_->sig.emit(type);
    }
};

std::shared_ptr<mock_pool> create_mock_pool(asio::io_context& ctx, pool_params&& params)
{
    return std::make_shared<mock_pool>(ctx.get_executor(), std::move(params));
}

class fixture : public io_context_fixture
{
private:
    // Pool (must be created using dynamic memory)
    std::shared_ptr<mock_pool> pool_;
    bool run_finished_{false};

public:
    fixture(pool_params&& params) : pool_(create_mock_pool(ctx, std::move(params)))
    {
        pool_->async_run([this](error_code ec) {
            run_finished_ = true;
            BOOST_TEST(ec == error_code());
        });
    }

    ~fixture()
    {
        // Finish the pool
        pool_->cancel();
        poll_until(ctx, &run_finished_);
    }

    mock_pool& pool() { return *pool_; }

    get_connection_task create_task(diagnostics* diag = nullptr, bool bind_slot = true)
    {
        return get_connection_task(*pool_, diag, bind_slot);
    }

    void check_shared_st(
        const diagnostics& expected_diag,
        std::size_t expected_num_pending,
        std::size_t expected_num_idle
    )
    {
        const auto& st = pool_->shared_state();
        BOOST_TEST(st.last_connect_diag == expected_diag);
        BOOST_TEST(st.num_pending_connections == expected_num_pending);
        BOOST_TEST(st.idle_list.size() == expected_num_idle);
    }

    // Waits for a status on a certain node
    void wait_for_status(
        mock_node& node,
        connection_status status,
        boost::source_location loc = BOOST_MYSQL_CURRENT_LOCATION
    )
    {
        poll_until(ctx, [&node, status]() { return node.status() == status; }, loc);
    }

    // Waits until there is at least num_nodes connections in the list
    void wait_for_num_nodes(std::size_t num_nodes, boost::source_location loc = BOOST_MYSQL_CURRENT_LOCATION)
    {
        poll_until(ctx, [this, num_nodes]() { return pool_->nodes().size() == num_nodes; }, loc);
    }

    // Wrapper for calling mock_connection::step()
    void step(
        mock_node& node,
        fn_type next_act,
        error_code ec = {},
        diagnostics diag = {},
        boost::source_location loc = BOOST_MYSQL_CURRENT_LOCATION
    )
    {
        node.connection().step(next_act, as_netresult, ec, diag).validate_no_error_nodiag(loc);
    }
};

// connection lifecycle
BOOST_AUTO_TEST_CASE(lifecycle_connect_error)
{
    // Setup
    const auto connect_diag = create_server_diag("Connection error!");
    const auto expected_diag = create_server_diag(
        "Last connection attempt failed with: er_aborting_connection [mysql.common-server:1152]: Connection "
        "error!"
    );
    pool_params params;
    params.retry_interval = std::chrono::seconds(2);
    fixture fix(std::move(params));

    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Connection trying to connect
    fix.wait_for_status(node, connection_status::connect_in_progress);
    fix.check_shared_st(diagnostics(), 1, 0);

    // Connect fails, so the connection goes to sleep. Diagnostics are stored in shared state.
    fix.step(node, fn_type::connect, common_server_errc::er_aborting_connection, connect_diag);
    fix.wait_for_status(node, connection_status::sleep_connect_failed_in_progress);
    fix.check_shared_st(expected_diag, 1, 0);

    // Advance until it's time to retry again
    mock_clock::advance_time_by(std::chrono::seconds(2));
    fix.wait_for_status(node, connection_status::connect_in_progress);
    fix.check_shared_st(expected_diag, 1, 0);

    // Connection connects successfully this time. Diagnostics have
    // been cleared and the connection is marked as idle
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);
    fix.check_shared_st(diagnostics(), 0, 1);
}

BOOST_AUTO_TEST_CASE(lifecycle_connect_timeout)
{
    // Setup
    pool_params params;
    params.connect_timeout = std::chrono::seconds(5);
    params.retry_interval = std::chrono::seconds(2);
    fixture fix(std::move(params));

    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Connection trying to connect
    fix.wait_for_status(node, connection_status::connect_in_progress);

    // Timeout elapses. Connect is considered failed
    const auto expected_diag = create_client_diag("Last connection attempt timed out");
    mock_clock::advance_time_by(std::chrono::seconds(5));
    fix.wait_for_status(node, connection_status::sleep_connect_failed_in_progress);
    fix.check_shared_st(expected_diag, 1, 0);

    // Advance until it's time to retry again
    mock_clock::advance_time_by(std::chrono::seconds(2));
    fix.wait_for_status(node, connection_status::connect_in_progress);
    fix.check_shared_st(expected_diag, 1, 0);

    // Connection connects successfully this time
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);
    fix.check_shared_st(diagnostics(), 0, 1);
}

BOOST_AUTO_TEST_CASE(lifecycle_return_without_reset)
{
    // Setup
    fixture fix(pool_params{});

    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Wait until a connection is successfully connected
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);
    fix.check_shared_st(diagnostics(), 0, 1);

    // Simulate a user picking the connection
    node.mark_as_in_use();
    fix.check_shared_st(diagnostics(), 0, 0);

    // Simulate a user returning the connection (without reset)
    fix.pool().return_connection(node, false);

    // The connection goes back to idle without invoking resets
    fix.wait_for_status(node, connection_status::idle);
    fix.check_shared_st(diagnostics(), 0, 1);
}

BOOST_AUTO_TEST_CASE(lifecycle_reset_success)
{
    // Setup
    fixture fix(pool_params{});

    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Wait until a connection is successfully connected, then pick it up
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);
    node.mark_as_in_use();

    // Simulate a user returning the connection (with reset)
    fix.pool().return_connection(node, true);

    // A reset is issued
    fix.wait_for_status(node, connection_status::reset_in_progress);
    fix.check_shared_st(diagnostics(), 1, 0);

    // Successful reset makes the connection idle again
    fix.step(node, fn_type::pipeline);
    fix.wait_for_status(node, connection_status::idle);
    fix.check_shared_st(diagnostics(), 0, 1);
}

BOOST_AUTO_TEST_CASE(lifecycle_reset_error)
{
    // Setup
    fixture fix(pool_params{});

    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Connect, pick up and return a connection
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);
    node.mark_as_in_use();
    fix.pool().return_connection(node, true);
    fix.wait_for_status(node, connection_status::reset_in_progress);

    // Reset fails. This triggers a reconnection. Diagnostics are not saved
    fix.step(node, fn_type::pipeline, common_server_errc::er_aborting_connection);
    fix.wait_for_status(node, connection_status::connect_in_progress);
    fix.check_shared_st(diagnostics(), 1, 0);

    // Reconnect succeeds. We're idle again
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);
    fix.check_shared_st(diagnostics(), 0, 1);
}

BOOST_AUTO_TEST_CASE(lifecycle_reset_timeout)
{
    // Setup
    pool_params params;
    params.ping_timeout = std::chrono::seconds(1);
    fixture fix(std::move(params));

    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Connect, pick up and return a connection
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);
    node.mark_as_in_use();
    fix.pool().return_connection(node, true);
    fix.wait_for_status(node, connection_status::reset_in_progress);

    // Reset times out. This triggers a reconnection
    mock_clock::advance_time_by(std::chrono::seconds(1));
    fix.wait_for_status(node, connection_status::connect_in_progress);
    fix.check_shared_st(diagnostics(), 1, 0);

    // Reconnect succeeds. We're idle again
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);
    fix.check_shared_st(diagnostics(), 0, 1);
}

BOOST_AUTO_TEST_CASE(lifecycle_reset_timeout_disabled)
{
    // Setup
    pool_params params;
    params.ping_timeout = std::chrono::seconds(0);
    fixture fix(std::move(params));

    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Connect, pick up and return a connection
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);
    node.mark_as_in_use();
    fix.pool().return_connection(node, true);
    fix.wait_for_status(node, connection_status::reset_in_progress);

    // Reset doesn't time out, regardless of how much time we wait
    mock_clock::advance_time_by(std::chrono::hours(9999));
    poll_until(fix.ctx, [&]() { return node.status() == connection_status::reset_in_progress; });
    fix.check_shared_st(diagnostics(), 1, 0);

    // Reset succeeds
    fix.step(node, fn_type::pipeline);
    fix.wait_for_status(node, connection_status::idle);
    fix.check_shared_st(diagnostics(), 0, 1);
}

BOOST_AUTO_TEST_CASE(lifecycle_ping_success)
{
    // Setup
    pool_params params;
    params.ping_interval = std::chrono::seconds(100);
    fixture fix(std::move(params));

    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Wait until a connection is successfully connected
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);

    // Wait until ping interval ellapses. This triggers a ping
    mock_clock::advance_time_by(std::chrono::seconds(100));
    fix.wait_for_status(node, connection_status::ping_in_progress);
    fix.check_shared_st(diagnostics(), 1, 0);

    // After ping succeeds, connection goes back to idle
    fix.step(node, fn_type::ping);
    fix.wait_for_status(node, connection_status::idle);
    fix.check_shared_st(diagnostics(), 0, 1);
}

BOOST_AUTO_TEST_CASE(lifecycle_ping_error)
{
    // Setup
    pool_params params;
    params.ping_interval = std::chrono::seconds(100);
    fixture fix(std::move(params));

    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Wait until a connection is successfully connected
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);

    // Wait until ping interval ellapses
    mock_clock::advance_time_by(std::chrono::seconds(100));

    // Ping fails. This triggers a reconnection. Diagnostics are not saved
    fix.step(node, fn_type::ping, common_server_errc::er_aborting_connection);
    fix.wait_for_status(node, connection_status::connect_in_progress);
    fix.check_shared_st(diagnostics(), 1, 0);

    // Reconnection succeeds
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);
    fix.check_shared_st(diagnostics(), 0, 1);
}

BOOST_AUTO_TEST_CASE(lifecycle_ping_timeout)
{
    // Setup
    pool_params params;
    params.ping_interval = std::chrono::seconds(100);
    params.ping_timeout = std::chrono::seconds(2);
    fixture fix(std::move(params));

    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Wait until a connection is successfully connected
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);

    // Wait until ping interval ellapses
    mock_clock::advance_time_by(std::chrono::seconds(100));
    fix.wait_for_status(node, connection_status::ping_in_progress);

    // Ping times out. This triggers a reconnection. Diagnostics are not saved
    mock_clock::advance_time_by(std::chrono::seconds(2));
    fix.wait_for_status(node, connection_status::connect_in_progress);

    // Reconnection succeeds
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);
    fix.check_shared_st(diagnostics(), 0, 1);
}

BOOST_AUTO_TEST_CASE(lifecycle_ping_timeout_disabled)
{
    // Setup
    pool_params params;
    params.ping_interval = std::chrono::seconds(100);
    params.ping_timeout = std::chrono::seconds(0);
    fixture fix(std::move(params));

    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Wait until a connection is successfully connected
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);

    // Wait until ping interval ellapses
    mock_clock::advance_time_by(std::chrono::seconds(100));
    fix.wait_for_status(node, connection_status::ping_in_progress);

    // Ping doesn't time out, regardless of how much we wait
    mock_clock::advance_time_by(std::chrono::hours(9999));
    poll_until(fix.ctx, [&]() { return node.status() == connection_status::ping_in_progress; });

    // Ping succeeds
    fix.step(node, fn_type::ping);
    fix.wait_for_status(node, connection_status::idle);
    fix.check_shared_st(diagnostics(), 0, 1);
}

BOOST_AUTO_TEST_CASE(lifecycle_ping_disabled)
{
    // Setup
    pool_params params;
    params.ping_interval = std::chrono::seconds(0);
    fixture fix(std::move(params));

    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Wait until a connection is successfully connected
    fix.step(node, fn_type::connect);
    fix.wait_for_status(node, connection_status::idle);

    // Connection won't ping, regardless of how much time we wait
    mock_clock::advance_time_by(std::chrono::hours(9999));
    poll_until(fix.ctx, [&]() { return node.status() == connection_status::idle; });
    fix.check_shared_st(diagnostics(), 0, 1);
}

// async_get_connection. Some cases are run with/without binding a cancel slot,
// and in thread-safe/unsafe mode
BOOST_AUTO_TEST_CASE(get_connection_immediate_completion)
{
    struct
    {
        bool bind_slot;
        bool thread_safe;
    } test_cases[] = {
        {false, false},
        {true,  false},
        {false, true },
        {true,  true },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT("bind_slot: " << tc.bind_slot << ", thread_safe: " << tc.thread_safe)
        {
            // Setup
            pool_params params;
            params.thread_safe = tc.thread_safe;
            fixture fix(std::move(params));

            fix.wait_for_num_nodes(1);
            auto& node = fix.pool().nodes().front();

            // Wait for a connection to be ready
            fix.step(node, fn_type::connect);
            fix.wait_for_status(node, connection_status::idle);

            // A request for a connection is issued. The request completes immediately.
            // In thread-safe mode, we still need to exit the strand, so we won't see an inline completion
            fix.create_task(nullptr, tc.bind_slot).wait(node, !tc.thread_safe);
            BOOST_TEST(node.status() == connection_status::in_use);
            BOOST_TEST(fix.pool().nodes().size() == 1u);
        }
    }
}

BOOST_AUTO_TEST_CASE(get_connection_wait_success)
{
    struct
    {
        bool bind_slot;
        bool thread_safe;
    } test_cases[] = {
        {false, false},
        {true,  false},
        {false, true },
        {true,  true },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT("bind_slot: " << tc.bind_slot << ", thread_safe: " << tc.thread_safe)
        {
            // Setup
            pool_params params;
            params.retry_interval = std::chrono::seconds(2);
            params.thread_safe = tc.thread_safe;
            fixture fix(std::move(params));

            fix.wait_for_num_nodes(1);
            auto& node = fix.pool().nodes().front();

            // Connection tries to connect and fails
            fix.step(node, fn_type::connect, common_server_errc::er_aborting_connection);
            fix.wait_for_status(node, connection_status::sleep_connect_failed_in_progress);

            // A request for a connection is issued. The request doesn't find
            // any available connection, and the current one is pending, so no new connections are created
            auto task = fix.create_task(nullptr, tc.bind_slot);
            fix.ctx.poll();
            BOOST_TEST(fix.pool().nodes().size() == 1u);

            // Retry interval ellapses and connection retries and succeeds
            mock_clock::advance_time_by(std::chrono::seconds(2));
            fix.step(node, fn_type::connect);

            // Request is fulfilled
            task.wait(node, false);
            BOOST_TEST(node.status() == connection_status::in_use);
            BOOST_TEST(fix.pool().nodes().size() == 1u);
        }
    }
}

BOOST_AUTO_TEST_CASE(get_connection_wait_op_cancelled)
{
    // verify that cancelling works for all cancel types and thread-safety modes
    struct
    {
        string_view name;
        asio::cancellation_type_t cancel_type;
        bool thread_safe;
    } test_cases[]{
        {"terminal",      asio::cancellation_type_t::terminal, false},
        {"partial",       asio::cancellation_type_t::partial,  false},
        {"total",         asio::cancellation_type_t::total,    false},
        {"safe_terminal", asio::cancellation_type_t::terminal, true },
        {"safe_partial",  asio::cancellation_type_t::partial,  true },
        {"safe_total",    asio::cancellation_type_t::total,    true },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            fixture fix(pool_params{});
            diagnostics diag;

            fix.wait_for_num_nodes(1);

            // A request for a connection is issued. The request doesn't find
            // any available connection, and the current one is pending, so no new connections are created
            auto task = fix.create_task(&diag);
            fix.ctx.poll();
            BOOST_TEST(fix.pool().nodes().size() == 1u);

            // The connection fails to connect
            fix.step(
                *fix.pool().nodes().begin(),
                fn_type::connect,
                common_server_errc::er_bad_db_error,
                create_server_diag("Bad db")
            );

            // The request gets cancelled. Appropriate diagnostics are returned
            task.cancel(tc.cancel_type);
            task.wait(client_errc::no_connection_available, false);
            BOOST_TEST(
                diag ==
                create_server_diag(
                    "Last connection attempt failed with: er_bad_db_error [mysql.common-server:1049]: Bad db"
                )
            );
            BOOST_TEST(fix.pool().nodes().size() == 1u);
        }
    }
}

BOOST_AUTO_TEST_CASE(get_connection_wait_op_cancelled_no_diag_available)
{
    // No diagnostics are returned if the first connect op hasn't finished yet
    // Setup
    fixture fix(pool_params{});
    diagnostics diag;

    // A request for a connection is issued
    auto task = fix.create_task(&diag);
    fix.ctx.poll();

    // The request gets cancelled. No diagnostics is available, so nothing is returned
    task.cancel();
    task.wait(client_errc::no_connection_available, false);
    BOOST_TEST(diag == diagnostics());
}

BOOST_AUTO_TEST_CASE(get_connection_wait_op_cancelled_timeout)
{
    // When connect times out, appropriate diagnostics are generated
    // Setup
    pool_params params;
    params.connect_timeout = std::chrono::seconds(5);
    fixture fix(std::move(params));
    diagnostics diag;

    fix.wait_for_num_nodes(1);

    // A request for a connection is issued. The request doesn't find
    // any available connection, and the current one is pending, so no new connections are created
    auto task = fix.create_task(&diag);
    fix.ctx.poll();
    BOOST_TEST(fix.pool().nodes().size() == 1u);

    // The connection attempt times out
    mock_clock::advance_time_by(std::chrono::seconds(6));
    fix.wait_for_status(fix.pool().nodes().front(), connection_status::sleep_connect_failed_in_progress);

    // The request gets cancelled. We get the expected error
    task.cancel();
    task.wait(client_errc::no_connection_available, false);
    BOOST_TEST(diag == create_client_diag("Last connection attempt timed out"));
    BOOST_TEST(fix.pool().nodes().size() == 1u);
}

BOOST_DATA_TEST_CASE_F(
    io_context_fixture,
    get_connection_wait_op_cancelled_pool_not_running,
    data::make({false, true})
)
{
    // If the op is cancelled because the pool is not running, appropriate diagnostics are issued
    // Setup
    pool_params params;
    params.thread_safe = sample;
    auto pool = create_mock_pool(ctx, std::move(params));
    diagnostics diag;

    // A request for a connection is issued. The request finds the pool not running,
    // and waits
    get_connection_task task(*pool, &diag);
    ctx.poll();

    // The request gets cancelled. We get the expected error
    task.cancel();
    task.wait(client_errc::pool_not_running, false);
    BOOST_TEST(diag == diagnostics());

    pool->cancel();
}

BOOST_AUTO_TEST_CASE(get_connection_wait_op_cancelled_diag_nullptr)
{
    // We don't crash if diag is nullptr
    // Setup
    fixture fix(pool_params{});

    fix.wait_for_num_nodes(1);

    // A request for a connection is issued. The request doesn't find
    // any available connection, and the current one is pending, so no new connections are created
    auto task = fix.create_task(nullptr);
    fix.ctx.poll();
    BOOST_TEST(fix.pool().nodes().size() == 1u);

    // The connection fails to connect
    fix.step(
        *fix.pool().nodes().begin(),
        fn_type::connect,
        common_server_errc::er_bad_db_error,
        create_server_diag("Bad db")
    );

    // The request gets cancelled
    task.cancel();
    task.wait(client_errc::no_connection_available, false);
    BOOST_TEST(fix.pool().nodes().size() == 1u);
}

BOOST_DATA_TEST_CASE(get_connection_wait_pool_cancelled, data::make({false, true}))
{
    // Setup
    pool_params params;
    params.thread_safe = sample;
    fixture fix(std::move(params));
    diagnostics diag;

    fix.wait_for_num_nodes(1);

    // A request for a connection is issued. The request doesn't find
    // any available connection, and the current one is pending, so no new connections are created
    auto task = fix.create_task(&diag);
    fix.ctx.poll();
    BOOST_TEST(fix.pool().nodes().size() == 1u);

    // The connection fails to connect
    fix.step(
        *fix.pool().nodes().begin(),
        fn_type::connect,
        common_server_errc::er_bad_db_error,
        create_server_diag("Bad db")
    );

    // The pool gets cancelled
    fix.pool().cancel();

    // This causes the request to get cancelled.
    // No diagnostics are provided here, as they're usually misleading
    task.wait(client_errc::pool_cancelled, false);
    BOOST_TEST(diag == diagnostics());
}

BOOST_AUTO_TEST_CASE(get_connection_connection_creation)
{
    // Setup
    pool_params params;
    params.initial_size = 1;
    params.max_size = 2;
    fixture fix(std::move(params));

    fix.wait_for_num_nodes(1);
    auto& node1 = fix.pool().nodes().front();

    // Wait for a connection to be ready, then get it from the pool
    fix.step(node1, fn_type::connect);
    fix.wait_for_status(node1, connection_status::idle);
    fix.create_task().wait(node1, true);

    // Another request is issued. The connection we have is in use, so another one is created.
    // Since this is not immediate, the task will need to wait
    auto task2 = fix.create_task();
    fix.ctx.poll();
    auto node2 = &*std::next(fix.pool().nodes().begin());

    // Connection connects successfully and is handed to us
    fix.step(*node2, fn_type::connect);
    task2.wait(*node2, false);
    BOOST_TEST(node2->status() == connection_status::in_use);
    BOOST_TEST(fix.pool().nodes().size() == 2u);

    // Another request is issued. All connections are in use but max size is already
    // reached, so no new connection is created
    auto task3 = fix.create_task();
    fix.ctx.poll();
    BOOST_TEST(fix.pool().nodes().size() == 2u);

    // When one of the connections is returned, the request is fulfilled
    fix.pool().return_connection(*node2, false);
    task3.wait(*node2, false);
    BOOST_TEST(fix.pool().nodes().size() == 2u);
}

BOOST_AUTO_TEST_CASE(get_connection_multiple_requests)
{
    // Setup
    pool_params params;
    params.initial_size = 2;
    params.max_size = 2;
    fixture fix(std::move(params));

    // 2 connection nodes are created from the beginning
    fix.wait_for_num_nodes(2);

    // Issue some parallel requests
    auto task1 = fix.create_task();
    auto task2 = fix.create_task();
    auto task3 = fix.create_task();
    auto task4 = fix.create_task(nullptr);
    auto task5 = fix.create_task();

    // Two connections can be created. These fulfill two requests
    auto node1 = &fix.pool().nodes().front();
    auto node2 = &*std::next(fix.pool().nodes().begin());
    fix.step(*node1, fn_type::connect);
    fix.step(*node2, fn_type::connect);
    task1.wait(*node1, false);
    task2.wait(*node2, false);

    // task4 gets cancelled
    task4.cancel();
    task4.wait(client_errc::no_connection_available, false);

    // A connection is returned. The first task to enter is served
    fix.pool().return_connection(*node1, true);
    fix.step(*node1, fn_type::pipeline);
    task3.wait(*node1, false);

    // The next connection to be returned is for task5
    fix.pool().return_connection(*node2, false);
    task5.wait(*node2, false);

    // Done
    BOOST_TEST(fix.pool().nodes().size() == 2u);
}

BOOST_AUTO_TEST_CASE(get_connection_supports_cancel_type)
{
    using ct = asio::cancellation_type_t;

    struct
    {
        string_view name;
        asio::cancellation_type_t input;
        bool expected;
    } test_cases[] = {
        {"none",                       ct::none,                               false},
        {"all",                        ct::all,                                true },
        {"terminal",                   ct::terminal,                           true },
        {"partial",                    ct::partial,                            true },
        {"total",                      ct::total,                              true },
        {"terminal | partial",         ct::terminal | ct::partial,             true },
        {"terminal | total",           ct::terminal | ct::total,               true },
        {"partial | total",            ct::partial | ct::total,                true },
        {"total | terminal | partial", ct::total | ct::terminal | ct::partial, true },
        {"unknown",                    static_cast<ct>(0xf000),                false},
        {"unknown | terminal",         static_cast<ct>(0xf000) | ct::terminal, true },
        {"unknown | total",            static_cast<ct>(0xf000) | ct::total,    true },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            BOOST_TEST(mock_pool::get_connection_supports_cancel_type(tc.input) == tc.expected);
        }
    }
}

// async_run
BOOST_AUTO_TEST_CASE(run_supports_cancel_type)
{
    using ct = asio::cancellation_type_t;

    struct
    {
        string_view name;
        asio::cancellation_type_t input;
        bool expected;
    } test_cases[] = {
        {"none",                       ct::none,                               false},
        {"all",                        ct::all,                                true },
        {"terminal",                   ct::terminal,                           true },
        {"partial",                    ct::partial,                            true },
        {"total",                      ct::total,                              false},
        {"terminal | partial",         ct::terminal | ct::partial,             true },
        {"terminal | total",           ct::terminal | ct::total,               true },
        {"partial | total",            ct::partial | ct::total,                true },
        {"total | terminal | partial", ct::total | ct::terminal | ct::partial, true },
        {"unknown",                    static_cast<ct>(0xf000),                false},
        {"unknown | terminal",         static_cast<ct>(0xf000) | ct::terminal, true },
        {"unknown | total",            static_cast<ct>(0xf000) | ct::total,    false},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            BOOST_TEST(mock_pool::run_supports_cancel_type(tc.input) == tc.expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(run_op_cancel)
{
    struct test_case
    {
        string_view name;
        asio::cancellation_type_t cancel_type;
        bool thread_safe;
    } test_cases[]{
        {"terminal",      asio::cancellation_type_t::terminal, false},
        {"partial",       asio::cancellation_type_t::partial,  false},
        {"safe_terminal", asio::cancellation_type_t::terminal, true },
        {"safe_partial",  asio::cancellation_type_t::partial,  true },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Create a pool
            io_context_fixture fix;
            pool_params params;
            params.thread_safe = tc.thread_safe;
            auto pool = create_mock_pool(fix.ctx, std::move(params));

            // Run with a bound signal
            bool run_finished = false;
            asio::cancellation_signal sig;
            pool->async_run(asio::bind_cancellation_slot(sig.slot(), [&run_finished](error_code ec) {
                run_finished = true;
                BOOST_TEST(ec == error_code());
            }));

            // Emit the signal. run should finish
            sig.emit(tc.cancel_type);
            poll_until(fix.ctx, &run_finished);

            // The pool has effectively been cancelled, as if cancel() had been called
            get_connection_task(*pool, nullptr).wait(client_errc::pool_cancelled, !tc.thread_safe);
        }
    }
}

// async_run with a bound signal that doesn't get emitted doesn't cause trouble
BOOST_FIXTURE_TEST_CASE(run_bound_signal, io_context_fixture)
{
    // Create a pool
    pool_params params;
    params.thread_safe = true;
    auto pool = create_mock_pool(ctx, std::move(params));

    // Run with a bound signal
    bool run_finished = false;
    asio::cancellation_signal sig;
    pool->async_run(asio::bind_cancellation_slot(sig.slot(), [&run_finished](error_code ec) {
        run_finished = true;
        BOOST_TEST(ec == error_code());
    }));

    // Cancel (but not using the signal)
    pool->cancel();

    // Finish successfully
    poll_until(ctx, &run_finished);
}

// pool size 0 works
BOOST_AUTO_TEST_CASE(initial_size_0)
{
    // Setup
    pool_params params;
    params.initial_size = 0;
    fixture fix(std::move(params));

    // No connections created at this point. A connection request arrives
    BOOST_TEST(fix.pool().nodes().size() == 0u);
    auto task = fix.create_task();

    // This creates a new connection, which fulfills the request
    fix.ctx.poll();
    BOOST_TEST(fix.pool().nodes().size() == 1u);
    fix.step(fix.pool().nodes().front(), fn_type::connect);
    task.wait(fix.pool().nodes().front(), false);
}

// pool_params have the intended effect
BOOST_AUTO_TEST_CASE(params_ssl_ctx_buffsize)
{
    // Setup. Pass a custom ssl context and buffer size
    // SSL context matching is performed using the underlying handle
    // because ssl::context provides no way to query the options previously set
    pool_params params;
    params.ssl_ctx.emplace(boost::asio::ssl::context::tlsv12_client);
    params.initial_buffer_size = 16u;
    auto handle = params.ssl_ctx->native_handle();
    fixture fix(std::move(params));

    // Wait for the node to be created
    fix.wait_for_num_nodes(1);

    // Check
    auto ctor_params = fix.pool().nodes().front().connection().ctor_params;
    BOOST_TEST_REQUIRE(ctor_params.ssl_context != nullptr);
    BOOST_TEST(ctor_params.ssl_context->native_handle() == handle);
    BOOST_TEST(ctor_params.initial_buffer_size == 16u);
}

BOOST_AUTO_TEST_CASE(params_connect_1)
{
    // Setup
    pool_params params;
    params.server_address.emplace_host_and_port("myhost", 1234);
    params.username = "myuser";
    params.password = "mypasswd";
    params.database = "mydb";
    params.ssl = boost::mysql::ssl_mode::disable;
    params.multi_queries = true;
    fixture fix(std::move(params));

    // Wait for the node to be created
    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Connect
    fix.step(node, fn_type::connect);

    // Check params
    const auto& cparams = node.connection().last_connect_params;
    BOOST_TEST(cparams.connection_collation == mysql_collations::utf8mb4_general_ci);
    BOOST_TEST(cparams.server_address.hostname() == "myhost");
    BOOST_TEST(cparams.server_address.port() == std::uint16_t(1234));
    BOOST_TEST(cparams.username == "myuser");
    BOOST_TEST(cparams.password == "mypasswd");
    BOOST_TEST(cparams.database == "mydb");
    BOOST_TEST(cparams.ssl == boost::mysql::ssl_mode::disable);
    BOOST_TEST(cparams.multi_queries == true);
}

BOOST_AUTO_TEST_CASE(params_connect_2)
{
    // Setup
    pool_params params;
    params.server_address.emplace_unix_path("/mysock");
    params.username = "myuser2";
    params.password = "mypasswd2";
    params.database = "mydb2";
    params.ssl = boost::mysql::ssl_mode::require;
    params.multi_queries = false;
    fixture fix(std::move(params));

    // Wait for the node to be created
    fix.wait_for_num_nodes(1);
    auto& node = fix.pool().nodes().front();

    // Connect
    fix.step(node, fn_type::connect);

    // Check params
    const auto& cparams = node.connection().last_connect_params;
    BOOST_TEST(cparams.connection_collation == mysql_collations::utf8mb4_general_ci);
    BOOST_TEST(cparams.server_address.unix_socket_path() == "/mysock");
    BOOST_TEST(cparams.username == "myuser2");
    BOOST_TEST(cparams.password == "mypasswd2");
    BOOST_TEST(cparams.database == "mydb2");
    BOOST_TEST(cparams.ssl == boost::mysql::ssl_mode::require);
    BOOST_TEST(cparams.multi_queries == false);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
