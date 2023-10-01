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
using boost::mysql::diagnostics;
using boost::mysql::error_code;

namespace {

struct sync_errc_maker
{
    static constexpr const char* name() { return "sync_errc"; }

    template <class Signature>
    struct type;

    template <class R, class Obj, class... Args>
    struct type<network_result<R>(Obj&, Args...)>
    {
        using impl = netfun_maker_sync_impl<R, Obj&, Args...>;
        using signature = std::function<network_result<R>(Obj&, Args...)>;
        using sync_sig = R (Obj::*)(Args..., error_code&, diagnostics&);

        static signature call(sync_sig sync) { return impl::sync_errc(sync); }
    };
};

}  // namespace

void boost::mysql::test::add_sync_errc(std::vector<er_network_variant*>& output)
{
    // Verify that all streams work
    static auto tcp = create_sync_variant<tcp_socket, sync_errc_maker>();
    static auto tcp_ssl = create_sync_variant<tcp_ssl_socket, sync_errc_maker>();
#if BOOST_ASIO_HAS_LOCAL_SOCKETS
    static auto unix = create_sync_variant<unix_socket, sync_errc_maker>();
    static auto unix_ssl = create_sync_variant<unix_ssl_socket, sync_errc_maker>();
#endif

    output.push_back(&tcp);
    output.push_back(&tcp_ssl);
#if BOOST_ASIO_HAS_LOCAL_SOCKETS
    output.push_back(&unix);
    output.push_back(&unix_ssl);
#endif
}
