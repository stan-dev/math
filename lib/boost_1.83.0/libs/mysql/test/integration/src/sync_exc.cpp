//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include "er_impl_common.hpp"
#include "test_common/netfun_helpers.hpp"
#include "test_integration/streams.hpp"

using namespace boost::mysql::test;

namespace {

struct sync_exc_maker
{
    static constexpr const char* name() { return "sync_exc"; }

    template <class Signature>
    struct type;

    template <class R, class Obj, class... Args>
    struct type<network_result<R>(Obj&, Args...)>
    {
        using impl = netfun_maker_sync_impl<R, Obj&, Args...>;
        using signature = std::function<network_result<R>(Obj&, Args...)>;
        using sync_sig = R (Obj::*)(Args...);

        static signature call(sync_sig sync) { return impl::sync_exc(sync); }
    };
};

}  // namespace

void boost::mysql::test::add_sync_exc(std::vector<er_network_variant*>& output)
{
    // Spotcheck
    static auto sync_exc = create_sync_variant<tcp_socket, sync_exc_maker>();
    output.push_back(&sync_exc);
}
