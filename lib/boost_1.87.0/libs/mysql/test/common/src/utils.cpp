//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/character_set.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/pipeline.hpp>
#include <boost/mysql/row.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/ssl_mode.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/access.hpp>

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/bind_executor.hpp>
#include <boost/asio/dispatch.hpp>
#include <boost/asio/execution/blocking.hpp>
#include <boost/asio/execution/relationship.hpp>
#include <boost/asio/execution_context.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/require.hpp>
#include <boost/assert/source_location.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <functional>
#include <iomanip>
#include <ostream>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/io_context_fixture.hpp"
#include "test_common/network_result.hpp"
#include "test_common/poll_until.hpp"
#include "test_common/printing.hpp"
#include "test_common/tracker_executor.hpp"
#include "test_common/validate_string_contains.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

//
// assert_buffer_equals.hpp
//
std::ostream& boost::mysql::test::operator<<(std::ostream& os, buffer_printer buff)
{
    os << std::setfill('0') << std::hex << "{ ";
    for (std::size_t i = 0; i < buff.buff.size(); ++i)
    {
        os << "0x" << std::setw(2) << static_cast<int>(buff.buff.data()[i]) << ", ";
    }
    return os << "}";
}

bool boost::mysql::test::buffer_equals(span<const std::uint8_t> b1, span<const std::uint8_t> b2)
{
    // If any of the buffers are empty (data() == nullptr), prevent
    // calling memcmp (UB)
    if (b1.size() == 0 || b2.size() == 0)
        return b1.size() == 0 && b2.size() == 0;

    if (b1.size() != b2.size())
        return false;

    return ::std::memcmp(b1.data(), b2.data(), b1.size()) == 0;
}

//
// printing.hpp
//

std::ostream& boost::mysql::operator<<(std::ostream& os, client_errc v) { return os << error_code(v); }

std::ostream& boost::mysql::operator<<(std::ostream& os, common_server_errc v) { return os << error_code(v); }

std::ostream& boost::mysql::operator<<(std::ostream& os, const diagnostics& diag)
{
    const auto& impl = detail::access::get_impl(diag);
    return os << "diagnostics{ " << (impl.is_server ? ".server_message" : ".client_message") << " = \""
              << impl.msg << "\" }";
}

std::ostream& boost::mysql::operator<<(std::ostream& os, const row_view& value)
{
    os << '{';
    if (!value.empty())
    {
        os << value[0];
        for (auto it = std::next(value.begin()); it != value.end(); ++it)
        {
            os << ", " << *it;
        }
    }
    return os << '}';
}

std::ostream& boost::mysql::operator<<(std::ostream& os, const row& r) { return os << row_view(r); }

static const char* to_string(metadata_mode v)
{
    switch (v)
    {
    case metadata_mode::full: return "metadata_mode::full";
    case metadata_mode::minimal: return "metadata_mode::minimal";
    default: return "<unknown metadata_mode>";
    }
}

std::ostream& boost::mysql::operator<<(std::ostream& os, metadata_mode v) { return os << ::to_string(v); }

static const char* to_string(ssl_mode v)
{
    switch (v)
    {
    case ssl_mode::disable: return "ssl_mode::disable";
    case ssl_mode::enable: return "ssl_mode::enable";
    case ssl_mode::require: return "ssl_mode::require";
    default: return "<unknown ssl_mode>";
    }
}

std::ostream& boost::mysql::operator<<(std::ostream& os, ssl_mode v) { return os << ::to_string(v); }

// character set
bool boost::mysql::operator==(const character_set& lhs, const character_set& rhs)
{
    if (lhs.name == nullptr || rhs.name == nullptr)
        return lhs.name == rhs.name;
    return std::strcmp(lhs.name, rhs.name) == 0 && lhs.next_char == rhs.next_char;
}

std::ostream& boost::mysql::operator<<(std::ostream& os, const character_set& v)
{
    if (v.name == nullptr)
        return os << "character_set()";
    else
        return os << "character_set(\"" << v.name << "\", .next_char? = " << static_cast<bool>(v.next_char)
                  << ")";
}

//
// tracker_executor.hpp
//

namespace {

// Are we in the call stack of an initiating function?
thread_local bool g_is_running_initiation = false;

// Produce unique executor IDs (start in 1)
std::atomic_int next_executor_id{1};

// The executor call stack
thread_local std::vector<int> g_executor_call_stack;

// Guard to remove an entry from the stack
struct executor_call_stack_guard
{
    executor_call_stack_guard(int executor_id) { g_executor_call_stack.push_back(executor_id); }
    executor_call_stack_guard(const executor_call_stack_guard&) = delete;
    executor_call_stack_guard(executor_call_stack_guard&&) = delete;
    executor_call_stack_guard& operator=(const executor_call_stack_guard&) = delete;
    executor_call_stack_guard& operator=(executor_call_stack_guard&&) = delete;
    ~executor_call_stack_guard() { g_executor_call_stack.pop_back(); }
};

// Wraps a function object to insert tracking of the caller executor
template <class Function>
struct tracker_executor_function
{
    int executor_id;
    Function fn;

    void operator()()
    {
        executor_call_stack_guard guard(executor_id);
        std::move(fn)();
    }
};

template <class Function>
tracker_executor_function<typename std::decay<Function>::type> create_tracker_executor_function(
    int executor_id,
    Function&& fn
)
{
    return {executor_id, std::forward<Function>(fn)};
}

}  // namespace

namespace boost {
namespace mysql {
namespace test {

class tracker_executor
{
public:
    tracker_executor(int id, boost::asio::any_io_executor ex) noexcept : id_(id), ex_(std::move(ex)) {}

    bool operator==(const tracker_executor& rhs) const noexcept { return id_ == rhs.id_ && ex_ == rhs.ex_; }
    bool operator!=(const tracker_executor& rhs) const noexcept { return !(*this == rhs); }
    int id() const { return id_; }

#ifdef BOOST_ASIO_USE_TS_EXECUTOR_AS_DEFAULT

    // TS-executor interface
    asio::execution_context& context() const noexcept { return ex_.context(); }
    void on_work_started() const noexcept { ex_.on_work_started(); }
    void on_work_finished() const noexcept { ex_.on_work_finished(); }

    template <typename Function, typename Allocator>
    void dispatch(Function&& f, const Allocator& a) const
    {
        ex_.dispatch(create_tracker_executor_function(id_, std::forward<Function>(f)), a);
    }

    template <typename Function, typename Allocator>
    void post(Function&& f, const Allocator& a) const
    {
        ex_.post(create_tracker_executor_function(id_, std::forward<Function>(f)), a);
    }

    template <typename Function, typename Allocator>
    void defer(Function&& f, const Allocator& a) const
    {
        ex_.defer(create_tracker_executor_function(id_, std::forward<Function>(f)), a);
    }

#else

    // Standard executors interface
    template <class Property>
    tracker_executor require(
        const Property& p,
        typename std::enable_if<asio::can_require<asio::any_io_executor, Property>::value>::type* = nullptr
    ) const
    {
        return tracker_executor(id_, asio::require(ex_, p));
    }

    template <class Property>
    tracker_executor prefer(
        const Property& p,
        typename std::enable_if<asio::can_prefer<asio::any_io_executor, Property>::value>::type* = nullptr
    ) const
    {
        return tracker_executor(id_, asio::prefer(ex_, p));
    }

    template <class Property>
    auto query(
        const Property& p,
        typename std::enable_if<asio::can_query<asio::any_io_executor, Property>::value>::type* = nullptr
    ) const -> decltype(asio::query(std::declval<boost::asio::any_io_executor>(), p))
    {
        return boost::asio::query(ex_, p);
    }

    template <typename Function>
    void execute(Function&& f) const
    {
        ex_.execute(create_tracker_executor_function(id_, std::forward<Function>(f)));
    }
#endif

private:
    int id_;
    boost::asio::any_io_executor ex_;
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

// Asio fails to detect that a type is equaility comparable under MSVC, so we need to do this
namespace boost {
namespace asio {
namespace traits {

template <>
struct equality_comparable<tracker_executor>
{
    static constexpr bool is_valid = true;
    static constexpr bool is_noexcept = true;
};

template <typename Function>
struct execute_member<tracker_executor, Function>
{
    static constexpr bool is_valid = true;
    static constexpr bool is_noexcept = false;
    using result_type = void;
};

template <class Property>
struct require_member<
    tracker_executor,
    Property,
    typename std::enable_if<asio::can_require<asio::any_io_executor, Property>::value>::type>
{
    static constexpr bool is_valid = true;
    static constexpr bool is_noexcept = false;
    using result_type = tracker_executor;
};

template <typename Property>
struct prefer_member<
    tracker_executor,
    Property,
    typename std::enable_if<asio::can_prefer<asio::any_io_executor, Property>::value>::type>
{
    static constexpr bool is_valid = true;
    static constexpr bool is_noexcept = false;
    using result_type = tracker_executor;
};

template <typename Property>
struct query_member<
    tracker_executor,
    Property,
    typename std::enable_if<asio::can_query<asio::any_io_executor, Property>::value>::type>
{
    static constexpr bool is_valid = true;
    static constexpr bool is_noexcept = false;
    using result_type = decltype(boost::asio::query(
        std::declval<const any_io_executor&>(),
        std::declval<const Property&>()
    ));
};

}  // namespace traits
}  // namespace asio
}  // namespace boost

boost::mysql::test::tracker_executor_result boost::mysql::test::create_tracker_executor(
    asio::any_io_executor inner
)
{
    int id = ++next_executor_id;
    return {id, tracker_executor(id, std::move(inner))};
}

boost::span<const int> boost::mysql::test::executor_stack() { return g_executor_call_stack; }

int boost::mysql::test::get_executor_id(asio::any_io_executor ex)
{
    auto* typed_ex = ex.target<tracker_executor>();
    return typed_ex ? typed_ex->id() : -1;
}
boost::mysql::test::initiation_guard::initiation_guard()
{
    BOOST_ASSERT(!g_is_running_initiation);
    g_is_running_initiation = true;
}

boost::mysql::test::initiation_guard::~initiation_guard()
{
    BOOST_ASSERT(g_is_running_initiation);
    g_is_running_initiation = false;
}

bool boost::mysql::test::is_initiation_function() { return g_is_running_initiation; }

//
// poll_until.hpp
//

void boost::mysql::test::poll_until(asio::io_context& ctx, const bool* done, source_location loc)
{
    poll_until(ctx, [done]() { return *done; }, loc);
}

void boost::mysql::test::poll_until(
    asio::io_context& ctx,
    const std::function<bool()>& done,
    source_location loc
)
{
    BOOST_TEST_CONTEXT("Called from " << loc)
    {
        using std::chrono::steady_clock;

        // Restart the context, in case it was stopped
        ctx.restart();

        // Poll until this time point
        constexpr std::chrono::seconds timeout(5);
        auto timeout_tp = steady_clock::now() + timeout;

        // Perform the polling
        while (!done() && steady_clock::now() < timeout_tp)
        {
            ctx.poll();
            std::this_thread::yield();
        }

        // Check for timeout
        BOOST_TEST_REQUIRE(done());
    }
}

void boost::mysql::test::run_in_context(
    asio::io_context& ctx,
    const std::function<void()>& fn,
    source_location loc
)
{
    bool finished = false;
    asio::dispatch(asio::bind_executor(ctx.get_executor(), [fn, &finished]() {
        fn();
        finished = true;
    }));
    poll_until(ctx, &finished, loc);
}

//
// io_context_fixture.hpp
//
boost::mysql::test::io_context_fixture::~io_context_fixture()
{
    // Verify that our tests don't leave unfinished work
    ctx.poll();
    BOOST_TEST(ctx.stopped());
}

//
// network_result.hpp
//

void boost::mysql::test::network_result_base::validate_immediate(bool expect_immediate, source_location loc)
    const
{
    BOOST_TEST_CONTEXT("Called from " << loc) { BOOST_TEST(was_immediate == expect_immediate); }
}

void boost::mysql::test::network_result_base::validate_no_error(source_location loc) const
{
    validate_error(error_code(), diagnostics(), loc);
}

void boost::mysql::test::network_result_base::validate_no_error_nodiag(source_location loc) const
{
    validate_error(error_code(), create_server_diag("<diagnostics unavailable>"), loc);
}

void boost::mysql::test::network_result_base::validate_error(
    error_code expected_err,
    const diagnostics& expected_diag,
    source_location loc
) const
{
    BOOST_TEST_CONTEXT("Called from " << loc)
    {
        BOOST_TEST(diag == expected_diag);
        BOOST_TEST_REQUIRE(err == expected_err);
    }
}

void boost::mysql::test::network_result_base::validate_error(
    common_server_errc expected_err,
    string_view expected_msg,
    source_location loc
)
{
    validate_error(expected_err, create_server_diag(expected_msg), loc);
}

void boost::mysql::test::network_result_base::validate_error(
    client_errc expected_err,
    string_view expected_msg,
    source_location loc
)
{
    validate_error(expected_err, create_client_diag(expected_msg), loc);
}

void boost::mysql::test::network_result_base::validate_error_contains(
    error_code expected_err,
    const std::vector<std::string>& pieces,
    source_location loc
)
{
    BOOST_TEST_CONTEXT("Called from " << loc)
    {
        validate_string_contains(diag.server_message(), pieces);
        BOOST_TEST_REQUIRE(err == expected_err);
    }
}

void boost::mysql::test::network_result_base::validate_any_error(source_location loc) const
{
    BOOST_TEST_CONTEXT("Called from " << loc) { BOOST_TEST_REQUIRE(err != error_code()); }
}

boost::mysql::test::test_detail::as_netres_handler_base::as_netres_handler_base(
    asio::io_context& ctx,
    asio::cancellation_slot slot,
    const diagnostics* diag
)
    : ex_(create_tracker_executor(ctx.get_executor())),
      immediate_ex_(create_tracker_executor(ctx.get_executor())),
      slot_(slot),
      diag_ptr(diag)
{
}

void boost::mysql::test::test_detail::as_netres_handler_base::complete_base(
    error_code ec,
    network_result_base& netres
) const
{
    // Are we in an immediate completion?
    bool is_immediate = is_initiation_function();

    // Check executor. The passed executor must be the top one in all cases.
    // Immediate completions must be dispatched through the immediate executor, too.
    // In all cases, we may encounter a bigger stack because of previous immediate completions.
    const std::array<int, 1> stack_data_regular{{ex_.executor_id}};
    const std::array<int, 2> stack_data_immediate{
        {immediate_ex_.executor_id, ex_.executor_id}
    };

    // Expected top of the executor stack
    boost::span<const int> expected_stack_top = is_immediate ? boost::span<const int>(stack_data_immediate)
                                                             : boost::span<const int>(stack_data_regular);

    // Actual top of the executor stack
    auto actual_stack_top = executor_stack().last(
        (std::min)(executor_stack().size(), expected_stack_top.size())
    );

    // Compare
    BOOST_TEST(actual_stack_top == expected_stack_top, boost::test_tools::per_element());

    // Assign error code and diagnostics
    netres.err = ec;
    if (diag_ptr)
        netres.diag = *diag_ptr;
    else
        netres.diag = create_server_diag("<diagnostics unavailable>");

    // Record immediate-ness
    netres.was_immediate = is_immediate;
}