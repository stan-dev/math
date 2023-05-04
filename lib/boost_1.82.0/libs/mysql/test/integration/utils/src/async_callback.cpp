//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/diagnostics.hpp>

#include <boost/asio/bind_executor.hpp>

#include "er_impl_common.hpp"
#include "netfun_helpers.hpp"
#include "streams.hpp"

using namespace boost::mysql::test;
using boost::mysql::diagnostics;

namespace {

struct async_callback_maker
{
    static constexpr const char* name() { return "async_callback"; }

    static void verify_tracked(const tracker_executor::tracked_values& v) { BOOST_TEST(v.total() > 0u); }

    template <class Signature>
    struct type;

    template <class R, class Obj, class... Args>
    struct type<network_result<R>(Obj&, Args...)>
    {
        using signature = std::function<network_result<R>(Obj&, Args...)>;
        using async_sig = void (Obj::*)(Args..., diagnostics&, bound_callback_token<R>&&);

        static signature call(async_sig fn)
        {
            return [fn](Obj& obj, Args... args) {
                tracker_executor::tracked_values tracked_vals;
                auto res = create_initial_netresult<R>();
                invoke_polyfill(
                    fn,
                    obj,
                    std::forward<Args>(args)...,
                    *res.diag,
                    boost::asio::bind_executor(
                        tracker_executor(get_context(obj.get_executor()), tracked_vals),
                        as_network_result<R>(res)
                    )
                );
                run_until_completion(obj.get_executor());
                verify_tracked(tracked_vals);
                return res;
            };
        }
    };
};

}  // namespace

void boost::mysql::test::add_async_callback(std::vector<er_network_variant*>& output)
{
    // Spotcheck for both streams
    static auto tcp = create_async_variant<tcp_socket, async_callback_maker>();
    static auto tcp_ssl = create_async_variant<tcp_ssl_socket, async_callback_maker>();
#if BOOST_ASIO_HAS_LOCAL_SOCKETS
    static auto unix_ssl = create_async_variant<unix_ssl_socket, async_callback_maker>();
#endif

    output.push_back(&tcp);
    output.push_back(&tcp_ssl);
#if BOOST_ASIO_HAS_LOCAL_SOCKETS
    output.push_back(&unix_ssl);
#endif
}
