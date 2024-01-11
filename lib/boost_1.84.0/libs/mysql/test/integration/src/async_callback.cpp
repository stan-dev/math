//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/connection.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/statement.hpp>

#include <boost/asio/bind_executor.hpp>

#include "er_impl_common.hpp"
#include "test_common/netfun_helpers.hpp"
#include "test_common/tracker_executor.hpp"
#include "test_integration/er_connection.hpp"
#include "test_integration/streams.hpp"

using namespace boost::mysql::test;

namespace boost {
namespace mysql {
namespace test {

template <class Stream>
class async_callback_connection : public connection_base<Stream>
{
    using conn_type = connection<Stream>;
    using base_type = connection_base<Stream>;

    template <class R, class... Args>
    using pmem_t = void (conn_type::*)(Args..., diagnostics&, as_network_result<R>&&);

    template <class R, class... Args>
    network_result<R> fn_impl(pmem_t<R, Args...> p, Args... args)
    {
        executor_info exec_info{};
        auto res = create_initial_netresult<R>();
        invoke_polyfill(
            p,
            this->conn(),
            std::forward<Args>(args)...,
            *res.diag,
            as_network_result<R>(res, create_tracker_executor(this->conn().get_executor(), &exec_info))
        );
        run_until_completion(this->conn().get_executor());
        BOOST_TEST(exec_info.total() > 0u);
        return res;
    }

public:
    BOOST_MYSQL_TEST_IMPLEMENT_ASYNC()
    static constexpr const char* name() noexcept { return "async_callback"; }
};

template <class Stream>
void add_async_callback_variant(std::vector<er_network_variant*>& output)
{
    add_variant<async_callback_connection<Stream>>(output);
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

void boost::mysql::test::add_async_callback(std::vector<er_network_variant*>& output)
{
    // Spotcheck for both streams
    add_async_callback_variant<tcp_socket>(output);
    add_async_callback_variant<tcp_ssl_socket>(output);
#if BOOST_ASIO_HAS_LOCAL_SOCKETS
    add_async_callback_variant<unix_socket>(output);
#endif
}
