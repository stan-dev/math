//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/asio/error.hpp>
#include <boost/asio/spawn.hpp>

#include <exception>

#include "er_impl_common.hpp"
#include "test_common/netfun_helpers.hpp"
#include "test_integration/streams.hpp"

using namespace boost::mysql::test;

namespace {

// Coroutines test async without diagnostics overloads
struct async_coroutine_maker
{
    static constexpr const char* name() { return "async_coroutines"; }

    template <class Signature>
    struct type;

    template <class R, class Obj, class... Args>
    struct type<network_result<R>(Obj&, Args...)>
    {
        using signature = std::function<network_result<R>(Obj&, Args...)>;
        using async_sig = R (Obj::*)(Args..., boost::asio::yield_context&&);

        static signature call(async_sig fn)
        {
            return [fn](Obj& obj, Args... args) {
                auto res = create_initial_netresult<R>(false);
                boost::asio::spawn(
                    obj.get_executor(),
                    [&](boost::asio::yield_context yield) {
                        invoke_and_assign(res, fn, obj, std::forward<Args>(args)..., yield[res.err]);
                    },
                    &rethrow_on_failure
                );
                run_until_completion(obj.get_executor());
                return res;
            };
        }
    };
};

}  // namespace

void boost::mysql::test::add_async_coroutines(std::vector<er_network_variant*>& output)
{
    static auto tcp = create_async_variant<tcp_socket, async_coroutine_maker>();
    output.push_back(&tcp);
}
