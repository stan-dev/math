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

namespace boost {
namespace mysql {
namespace test {

template <class Stream>
class sync_exc_connection : public connection_base<Stream>
{
    using conn_type = connection<Stream>;
    using base_type = connection_base<Stream>;

    template <class R, class... Args>
    using pmem_t = R (conn_type::*)(Args...);

    template <class R, class... Args>
    network_result<R> fn_impl(pmem_t<R, Args...> p, Args... args)
    {
        network_result<R> res;
        try
        {
            invoke_and_assign(res, p, this->conn(), std::forward<Args>(args)...);
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
    static constexpr const char* name() noexcept { return "sync_exc"; }
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

void boost::mysql::test::add_sync_exc(std::vector<er_network_variant*>& output)
{
    // Spotcheck
    add_variant<sync_exc_connection<tcp_socket>>(output);
}
