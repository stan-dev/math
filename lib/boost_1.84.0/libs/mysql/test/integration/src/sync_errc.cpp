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
using boost::mysql::error_code;

namespace boost {
namespace mysql {
namespace test {

template <class Stream>
class sync_errc_connection : public connection_base<Stream>
{
    using conn_type = connection<Stream>;
    using base_type = connection_base<Stream>;

    // workaround for gcc5
    template <class R, class... Args>
    struct pmem
    {
        using type = R (conn_type::*)(Args..., error_code&, diagnostics&);
    };

    template <class R, class... Args>
    network_result<R> fn_impl(typename pmem<R, Args...>::type p, Args... args)
    {
        auto res = create_initial_netresult<R>();
        invoke_and_assign(res, p, this->conn(), std::forward<Args>(args)..., res.err, *res.diag);
        return res;
    }

public:
// MSVC complains about passing empty tokens, which is valid C++
#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable : 4003)
#endif
    BOOST_MYSQL_TEST_IMPLEMENT_SYNC()
#ifdef BOOST_MSVC
#pragma warning(pop)
#endif
    static constexpr const char* name() noexcept { return "sync_errc"; }
};

template <class Stream>
void add_sync_errc_variant(std::vector<er_network_variant*>& output)
{
    add_variant<sync_errc_connection<Stream>>(output);
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

void boost::mysql::test::add_sync_errc(std::vector<er_network_variant*>& output)
{
    // Verify that all streams work
    add_sync_errc_variant<tcp_socket>(output);
    add_sync_errc_variant<tcp_ssl_socket>(output);
#if BOOST_ASIO_HAS_LOCAL_SOCKETS
    add_sync_errc_variant<unix_socket>(output);
    add_sync_errc_variant<unix_ssl_socket>(output);
#endif
}
