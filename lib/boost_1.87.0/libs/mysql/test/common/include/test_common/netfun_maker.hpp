//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_NETFUN_MAKER_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_NETFUN_MAKER_HPP

#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>

#include <boost/system/system_error.hpp>

#include <type_traits>

#include "test_common/create_diagnostics.hpp"
#include "test_common/network_result.hpp"

namespace boost {
namespace mysql {
namespace test {

template <class Fn, class... Args>
auto invoke_polyfill(Fn fn, Args&&... args) ->
    typename std::enable_if<
        std::is_function<typename std::remove_pointer<Fn>::type>::value,
        decltype(fn(std::forward<Args>(args)...))>::type
{
    return fn(std::forward<Args>(args)...);
}

template <class Pmem, class Obj, class... Args>
auto invoke_polyfill(Pmem fn, Obj& obj, Args&&... args) ->
    typename std::enable_if<
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

template <class R, class IOObject, class... Args>
struct netfun_maker_impl
{
    using signature = std::function<network_result<R>(IOObject&, Args...)>;

    template <class Pfn>
    static signature sync_errc(Pfn fn)
    {
        return [fn](IOObject& obj, Args... args) {
            network_result<R> res{
                common_server_errc::er_no,
                create_server_diag("diagnostics not cleared properly")
            };
            invoke_and_assign(res, fn, obj, std::forward<Args>(args)..., res.err, res.diag);
            return res;
        };
    }

    template <class Pfn>
    static signature sync_errc_nodiag(Pfn fn)
    {
        return [fn](IOObject& obj, Args... args) {
            network_result<R> res{common_server_errc::er_no, create_server_diag("<diagnostics unavailable>")};
            invoke_and_assign(res, fn, obj, std::forward<Args>(args)..., res.err);
            return res;
        };
    }

    template <class Pfn>
    static signature sync_exc(Pfn fn)
    {
        return [fn](IOObject& obj, Args... args) {
            network_result<R> res;
            try
            {
                invoke_and_assign(res, fn, obj, std::forward<Args>(args)...);
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

    template <class Pfn>
    static signature async_diag(Pfn fn)
    {
        return [fn](IOObject& obj, Args... args) {
            diagnostics diag;  // checks for clearing diag are performed by as_netresult
            return invoke_polyfill(fn, obj, std::forward<Args>(args)..., diag, as_netresult).run();
        };
    }

    template <class Pfn>
    static signature async_nodiag(Pfn fn)
    {
        return [fn](IOObject& obj, Args... args) {
            return invoke_polyfill(fn, obj, std::forward<Args>(args)..., as_netresult).run();
        };
    }
};

template <class R, class Obj, class... Args>
class netfun_maker
{
    using impl = netfun_maker_impl<R, Obj&, Args...>;

public:
    using signature = std::function<network_result<R>(Obj&, Args...)>;
    using sig_sync_errc = R (Obj::*)(Args..., error_code&, diagnostics&);
    using sig_sync_errc_nodiag = R (Obj::*)(Args..., error_code&);
    using sig_sync_errc_nodiag_old =
        error_code (Obj::*)(Args..., error_code&);  // support old Asio signatures
    using sig_sync_exc = R (Obj::*)(Args...);
    using sig_async_diag = runnable_network_result<R> (Obj::*)(Args..., diagnostics&, const as_netresult_t&);
    using sig_async_nodiag = runnable_network_result<R> (Obj::*)(Args..., const as_netresult_t&);

    static signature sync_errc(sig_sync_errc pfn) { return impl::sync_errc(pfn); }
    static signature sync_errc_nodiag(sig_sync_errc_nodiag pfn) { return impl::sync_errc_nodiag(pfn); }
    static signature sync_errc_nodiag(sig_sync_errc_nodiag_old pfn) { return impl::sync_errc_nodiag(pfn); }
    static signature sync_exc(sig_sync_exc pfn) { return impl::sync_exc(pfn); }
    static signature async_diag(sig_async_diag pfn) { return impl::async_diag(pfn); }
    static signature async_nodiag(sig_async_nodiag pfn) { return impl::async_nodiag(pfn); }
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
