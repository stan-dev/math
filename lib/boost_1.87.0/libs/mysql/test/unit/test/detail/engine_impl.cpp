//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>

#include <boost/mysql/detail/any_resumable_ref.hpp>
#include <boost/mysql/detail/engine_impl.hpp>
#include <boost/mysql/detail/next_action.hpp>

#include <boost/asio/any_completion_handler.hpp>
#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/async_result.hpp>
#include <boost/asio/bind_executor.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/coroutine.hpp>
#include <boost/asio/deferred.hpp>
#include <boost/asio/dispatch.hpp>
#include <boost/asio/error.hpp>
#include <boost/asio/post.hpp>
#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>

#include "test_common/create_diagnostics.hpp"
#include "test_common/io_context_fixture.hpp"
#include "test_common/netfun_maker.hpp"
#include "test_common/network_result.hpp"
#include "test_common/poll_until.hpp"
#include "test_unit/printing.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql::detail;
namespace asio = boost::asio;
using boost::mysql::error_code;

BOOST_AUTO_TEST_SUITE(test_engine_impl)

// Satisfies the EngineStream concept
class mock_engine_stream
{
    asio::any_io_executor ex_;
    error_code op_error_;  // Operations complete with this error code

    std::size_t size_or_zero(std::size_t sz) const { return op_error_ ? 0u : sz; }

    void record_read_call(asio::mutable_buffer buff, bool use_ssl)
    {
        calls.push_back(next_action::read({
            {static_cast<std::uint8_t*>(buff.data()), buff.size()},
            use_ssl
        }));
    }

    void record_write_call(asio::const_buffer buff, bool use_ssl)
    {
        calls.push_back(next_action::write({
            {static_cast<const std::uint8_t*>(buff.data()), buff.size()},
            use_ssl
        }));
    }

    template <class CompletionToken>
    void complete_immediate(CompletionToken&& token)
    {
        asio::post(ex_, asio::deferred([this]() {
                       return asio::deferred.values(op_error_);
                   }))(std::forward<CompletionToken>(token));
    }

    template <class CompletionToken>
    void complete_immediate(CompletionToken&& token, std::size_t bytes)
    {
        asio::post(ex_, asio::deferred([this, bytes]() {
                       return asio::deferred.values(op_error_, size_or_zero(bytes));
                   }))(std::forward<CompletionToken>(token));
    }

public:
    std::vector<next_action> calls;

    mock_engine_stream(asio::any_io_executor ex, error_code op_error = error_code())
        : ex_(std::move(ex)), op_error_(op_error)
    {
    }

    using executor_type = asio::any_io_executor;
    executor_type get_executor() { return ex_; }

    bool supports_ssl() const { return true; }

    void set_endpoint(const void*) {}

    // Reading
    std::size_t read_some(asio::mutable_buffer buff, bool use_ssl, error_code& ec)
    {
        record_read_call(buff, use_ssl);
        ec = op_error_;
        return size_or_zero(buff.size());
    }

    template <class CompletionToken>
    void async_read_some(asio::mutable_buffer buff, bool use_ssl, CompletionToken&& token)
    {
        record_read_call(buff, use_ssl);
        complete_immediate(std::forward<CompletionToken>(token), buff.size());
    }

    // Writing
    std::size_t write_some(asio::const_buffer buff, bool use_ssl, error_code& ec)
    {
        record_write_call(buff, use_ssl);
        ec = op_error_;
        return size_or_zero(buff.size());
    }

    template <class CompletionToken>
    void async_write_some(asio::const_buffer buff, bool use_ssl, CompletionToken&& token)
    {
        record_write_call(buff, use_ssl);
        complete_immediate(std::forward<CompletionToken>(token), buff.size());
    }

    // SSL
    void ssl_handshake(error_code& ec)
    {
        calls.push_back(next_action::ssl_handshake());
        ec = op_error_;
    }

    template <class CompletionToken>
    void async_ssl_handshake(CompletionToken&& token)
    {
        calls.push_back(next_action::ssl_handshake());
        complete_immediate(std::forward<CompletionToken>(token));
    }

    void ssl_shutdown(error_code& ec)
    {
        calls.push_back(next_action::ssl_shutdown());
        ec = op_error_;
    }

    template <class CompletionToken>
    void async_ssl_shutdown(CompletionToken&& token)
    {
        calls.push_back(next_action::ssl_shutdown());
        complete_immediate(std::forward<CompletionToken>(token));
    }

    // Connect and close
    void connect(error_code& ec)
    {
        calls.push_back(next_action::connect());
        ec = op_error_;
    }

    template <class CompletionToken>
    void async_connect(CompletionToken&& token)
    {
        calls.push_back(next_action::connect());
        complete_immediate(std::forward<CompletionToken>(token));
    }

    void close(error_code& ec)
    {
        calls.push_back(next_action::close());
        ec = op_error_;
    }
};

// Helpers to run the sync and async versions uniformly.
// engine_impl must be adapted because async_run() takes a completion handler,
// rather than a generic completion token
struct test_engine
{
    engine_impl<mock_engine_stream> value;

    struct initiation_t
    {
        asio::any_io_executor ex;

        // Required for compatibility with as_netresult
        using executor_type = asio::any_io_executor;
        executor_type get_executor() const { return ex; }

        void operator()(
            asio::any_completion_handler<void(error_code)> handler,
            engine* eng,
            any_resumable_ref resumable
        ) const
        {
            eng->async_run(resumable, std::move(handler));
        }
    };

    void run(any_resumable_ref resumable, error_code& err) { value.run(resumable, err); }

    template <class CompletionToken>
    auto async_run(any_resumable_ref resumable, CompletionToken&& token)
        -> decltype(asio::async_initiate<CompletionToken, void(error_code)>(
            initiation_t{value.get_executor()},
            token,
            &value,
            resumable
        ))
    {
        return asio::async_initiate<CompletionToken, void(error_code)>(
            initiation_t{value.get_executor()},
            token,
            &value,
            resumable
        );
    }
};

using netmaker_t = netfun_maker<void, test_engine, any_resumable_ref>;
const auto sync_fn = netmaker_t::sync_errc_nodiag(&test_engine::run);
const auto async_fn = netmaker_t::async_nodiag(&test_engine::async_run);
using signature_t = netmaker_t::signature;

// A mock for a sans-io algorithm. Can be converted to any_resumable_ref
struct mock_algo
{
    using call_args = std::pair<error_code, std::size_t>;

    std::size_t current_call{0};
    std::vector<next_action> acts;
    std::vector<call_args> calls;

    mock_algo(next_action act) : acts{act} {}
    mock_algo(next_action act1, next_action act2) : acts{act1, act2} {}

    void check_calls(const std::vector<call_args>& expected) { BOOST_TEST(calls == expected); }

    next_action resume(error_code ec, std::size_t bytes_transferred)
    {
        calls.push_back(call_args{ec, bytes_transferred});
        auto idx = current_call++;
        return idx < acts.size() ? acts[idx] : next_action();
    }
};

// returning next_action::read calls the relevant stream function
BOOST_AUTO_TEST_CASE(next_action_read)
{
    struct
    {
        const char* name;
        signature_t fn;
        bool ssl_active;
    } test_cases[] = {
        {"sync_ssl_active",    sync_fn,  true },
        {"sync_ssl_inactive",  sync_fn,  false},
        {"async_ssl_active",   async_fn, true },
        {"async_ssl_inactive", async_fn, false},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            io_context_fixture fix;
            std::array<std::uint8_t, 8> buff{};
            mock_algo algo(next_action::read({buff, tc.ssl_active}));
            test_engine eng{fix.ctx.get_executor()};

            tc.fn(eng, any_resumable_ref(algo)).validate_no_error_nodiag();
            BOOST_TEST(eng.value.stream().calls.size() == 1u);
            BOOST_TEST(eng.value.stream().calls[0].type() == next_action_type::read);
            BOOST_TEST(eng.value.stream().calls[0].read_args().use_ssl == tc.ssl_active);
            BOOST_TEST(eng.value.stream().calls[0].read_args().buffer.data() == buff.data());
            BOOST_TEST(eng.value.stream().calls[0].read_args().buffer.size() == buff.size());
            algo.check_calls({
                {error_code(), 0u},
                {error_code(), 8u}
            });
            // The testing infrastructure checks that we post correctly in async functions
        }
    }
}

// returning next_action::write calls the relevant stream function
BOOST_AUTO_TEST_CASE(next_action_write)
{
    struct
    {
        const char* name;
        signature_t fn;
        bool ssl_active;
    } test_cases[] = {
        {"sync_ssl_active",    sync_fn,  true },
        {"sync_ssl_inactive",  sync_fn,  false},
        {"async_ssl_active",   async_fn, true },
        {"async_ssl_inactive", async_fn, false},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            io_context_fixture fix;
            const std::array<std::uint8_t, 4> buff{};
            mock_algo algo(next_action::write({buff, tc.ssl_active}));
            test_engine eng{fix.ctx.get_executor()};

            tc.fn(eng, any_resumable_ref(algo)).validate_no_error_nodiag();
            BOOST_TEST(eng.value.stream().calls.size() == 1u);
            BOOST_TEST(eng.value.stream().calls[0].type() == next_action_type::write);
            BOOST_TEST(eng.value.stream().calls[0].write_args().use_ssl == tc.ssl_active);
            BOOST_TEST(eng.value.stream().calls[0].write_args().buffer.data() == buff.data());
            BOOST_TEST(eng.value.stream().calls[0].write_args().buffer.size() == buff.size());
            algo.check_calls({
                {error_code(), 0u},
                {error_code(), 4u}
            });
            // The testing infrastructure checks that we post correctly in async functions
        }
    }
}

// returning next_action::connect/ssl_handshake/ssl_shutdown/close calls the relevant stream function
BOOST_AUTO_TEST_CASE(next_action_other)
{
    struct
    {
        const char* name;
        signature_t fn;
        next_action act;
    } test_cases[] = {
        {"connect_sync",        sync_fn,  next_action::connect()      },
        {"connect_async",       async_fn, next_action::connect()      },
        {"ssl_handshake_sync",  sync_fn,  next_action::ssl_handshake()},
        {"ssl_handshake_async", async_fn, next_action::ssl_handshake()},
        {"ssl_shutdown_sync",   sync_fn,  next_action::ssl_shutdown() },
        {"ssl_shutdown_async",  async_fn, next_action::ssl_shutdown() },
        {"close_sync",          sync_fn,  next_action::close()        },
        {"close_async",         async_fn, next_action::close()        },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            io_context_fixture fix;
            mock_algo algo(tc.act);
            test_engine eng{fix.ctx.get_executor()};

            tc.fn(eng, any_resumable_ref(algo)).validate_no_error_nodiag();
            BOOST_TEST(eng.value.stream().calls.size() == 1u);
            BOOST_TEST(eng.value.stream().calls[0].type() == tc.act.type());
            algo.check_calls({
                {error_code(), 0u},
                {error_code(), 0u}
            });
            // The testing infrastructure checks that we post correctly in async functions
        }
    }
}

// Stream errors get propagated to the algorithm and don't exit the loop
BOOST_AUTO_TEST_CASE(stream_errors)
{
    std::array<std::uint8_t, 8> buff{};
    const std::array<std::uint8_t, 4> cbuff{};

    struct
    {
        const char* name;
        signature_t fn;
        next_action act;
    } test_cases[] = {
        {"read_sync",           sync_fn,  next_action::read({buff, false})  },
        {"read_async",          async_fn, next_action::read({buff, false})  },
        {"write_sync",          sync_fn,  next_action::write({cbuff, false})},
        {"write_async",         async_fn, next_action::write({cbuff, false})},
        {"connect_sync",        sync_fn,  next_action::connect()            },
        {"connect_async",       async_fn, next_action::connect()            },
        {"ssl_handshake_sync",  sync_fn,  next_action::ssl_handshake()      },
        {"ssl_handshake_async", async_fn, next_action::ssl_handshake()      },
        {"ssl_shutdown_sync",   sync_fn,  next_action::ssl_shutdown()       },
        {"ssl_shutdown_async",  async_fn, next_action::ssl_shutdown()       },
        {"close_sync",          sync_fn,  next_action::close()              },
        {"close_async",         async_fn, next_action::close()              },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            io_context_fixture fix;
            mock_algo algo(tc.act);
            test_engine eng{
                {fix.ctx.get_executor(), asio::error::already_open}
            };

            tc.fn(eng, any_resumable_ref(algo))
                .validate_no_error_nodiag();  // Error gets swallowed by the algo
            BOOST_TEST(eng.value.stream().calls.size() == 1u);
            BOOST_TEST(eng.value.stream().calls[0].type() == tc.act.type());
            algo.check_calls({
                {error_code(),              0u},
                {asio::error::already_open, 0u}
            });
            // The testing infrastructure checks that we post correctly in async functions
        }
    }
}

// Returning an error or next_action() from resume in the first call works correctly
BOOST_AUTO_TEST_CASE(resume_error_immediate_sync)
{
    struct
    {
        const char* name;
        error_code ec;
    } test_cases[] = {
        {"success", error_code()        },
        {"error",   asio::error::no_data},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            io_context_fixture fix;
            mock_algo algo(next_action(tc.ec));
            test_engine eng{fix.ctx.get_executor()};

            sync_fn(eng, any_resumable_ref(algo))
                .validate_error(tc.ec, create_server_diag("<diagnostics unavailable>"));
            BOOST_TEST(eng.value.stream().calls.size() == 0u);
            algo.check_calls({
                {error_code(), 0u}
            });
        }
    }
}

BOOST_AUTO_TEST_CASE(resume_error_immediate_async)
{
    struct
    {
        const char* name;
        error_code ec;
    } test_cases[] = {
        {"success", error_code()        },
        {"error",   asio::error::no_data},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // We need to call the initiation function from a context thread
            // to get immediate completions
            io_context_fixture fix;
            run_in_context(fix.ctx, [&]() {
                // Setup
                mock_algo algo(next_action(tc.ec));
                test_engine eng{fix.ctx.get_executor()};

                // Run the function
                eng.async_run(any_resumable_ref(algo), as_netresult)
                    .run()
                    .validate_immediate(true)
                    .validate_error(tc.ec, create_server_diag("<diagnostics unavailable>"));

                // Check
                BOOST_TEST(eng.value.stream().calls.size() == 0u);
                algo.check_calls({
                    {error_code(), 0u}
                });
            });
        }
    }
}

// Returning an error or next_action() from resume in successive calls works correctly
BOOST_AUTO_TEST_CASE(resume_error_successive_calls)
{
    struct
    {
        const char* name;
        signature_t fn;
        error_code ec;
    } test_cases[] = {
        {"success_sync",  sync_fn,  error_code()        },
        {"success_async", async_fn, error_code()        },
        {"error_sync",    sync_fn,  asio::error::no_data},
        {"error_async",   async_fn, asio::error::no_data},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            io_context_fixture fix;
            mock_algo algo(next_action::connect(), next_action(tc.ec));
            test_engine eng{fix.ctx.get_executor()};

            tc.fn(eng, any_resumable_ref(algo))
                .validate_error(tc.ec, create_server_diag("<diagnostics unavailable>"));
            BOOST_TEST(eng.value.stream().calls.size() == 1u);
            BOOST_TEST(eng.value.stream().calls[0].type() == next_action_type::connect);
            algo.check_calls({
                {error_code(), 0u},
                {error_code(), 0u}
            });
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
