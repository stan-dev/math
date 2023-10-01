//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_NETFUN_HELPERS_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_NETFUN_HELPERS_HPP

#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>

#include <boost/mysql/detail/void_t.hpp>

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/system/system_error.hpp>

#include <cstddef>
#include <functional>
#include <type_traits>

#include "test_common/create_diagnostics.hpp"
#include "test_common/network_result.hpp"

// Helper functions and classes to implement netmakers
// (the insfrastructure to run sync and async code as parameterized tests)
// for both unit and integ tests

namespace boost {
namespace mysql {
namespace test {

// Completion callback that saves result into a network_result
// and binds it to an executor
template <class R>
class as_network_result
{
    network_result<R>* netresult_;
    asio::any_io_executor executor_;

public:
    as_network_result(network_result<R>& netresult, asio::any_io_executor exec)
        : netresult_(&netresult), executor_(exec)
    {
    }

    using executor_type = asio::any_io_executor;
    asio::any_io_executor get_executor() const { return executor_; }

    void operator()(error_code ec) const noexcept { netresult_->err = ec; }

    template <class Arg>
    void operator()(error_code ec, Arg&& arg) const noexcept
    {
        netresult_->err = ec;
        netresult_->value = std::forward<Arg>(arg);
    }
};

template <class Fn, class... Args>
auto invoke_polyfill(Fn fn, Args&&... args) -> typename std::enable_if<
    std::is_function<typename std::remove_pointer<Fn>::type>::value,
    decltype(fn(std::forward<Args>(args)...))>::type
{
    return fn(std::forward<Args>(args)...);
}

template <class Pmem, class Obj, class... Args>
auto invoke_polyfill(Pmem fn, Obj& obj, Args&&... args) -> typename std::enable_if<
    std::is_member_function_pointer<Pmem>::value,
    decltype((obj.*fn)(std::forward<Args>(args)...))>::type
{
    return (obj.*fn)(std::forward<Args>(args)...);
}

template <class T, class... InvokeArgs>
void invoke_and_assign(network_result<T>& output, InvokeArgs&&... args)
{
    output.value = invoke_polyfill(std::forward<InvokeArgs>(args)...);
}

template <class... InvokeArgs>
void invoke_and_assign(network_result<void>&, InvokeArgs&&... args)
{
    invoke_polyfill(std::forward<InvokeArgs>(args)...);
}

template <class R>
network_result<R> create_initial_netresult(bool with_diag = true)
{
    network_result<R> res(boost::mysql::make_error_code(boost::mysql::common_server_errc::er_no));
    if (with_diag)
        res.diag = create_server_diag("diagnostics not cleared properly");
    return res;
}

// The synchronous implementations are common between unit and integ tests
template <class R, class... Args>
struct netfun_maker_sync_impl
{
    using signature = std::function<network_result<R>(Args...)>;

    template <class Pfn>
    static signature sync_errc(Pfn fn)
    {
        return [fn](Args... args) {
            auto res = create_initial_netresult<R>();
            invoke_and_assign(res, fn, std::forward<Args>(args)..., res.err, *res.diag);
            return res;
        };
    }

    template <class Pfn>
    static signature sync_exc(Pfn fn)
    {
        return [fn](Args... args) {
            network_result<R> res;
            try
            {
                invoke_and_assign(res, fn, std::forward<Args>(args)...);
            }
            catch (const boost::mysql::error_with_diagnostics& err)
            {
                res.err = err.code();
                res.diag = err.get_diagnostics();
            }
            catch (const boost::system::system_error& err)
            {
                res.err = err.code();
            }
            return res;
        };
    }
};

inline boost::asio::io_context& get_context(boost::asio::any_io_executor ex) noexcept
{
    return static_cast<boost::asio::io_context&>(ex.context());
}

inline void run_until_completion(boost::asio::any_io_executor ex)
{
    auto& ctx = get_context(ex);
    ctx.restart();
    ctx.run();
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
