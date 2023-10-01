//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_UNIT_NETFUN_MAKER_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_UNIT_NETFUN_MAKER_HPP

#include <boost/asio/any_io_executor.hpp>

#include "test_common/netfun_helpers.hpp"
#include "test_common/tracker_executor.hpp"

namespace boost {
namespace mysql {
namespace test {

// Note: the async implementations are only valid for unit tests.
// These rely on the stream dispatching completion handlers directly
// via post(), which only happens with test_stream.
template <class R, class IOObject, class... Args>
struct netfun_maker_impl : netfun_maker_sync_impl<R, IOObject, Args...>
{
    using signature = std::function<network_result<R>(IOObject&, Args...)>;

    template <class Pfn>
    static signature async_errinfo(Pfn fn)
    {
        return [fn](IOObject& obj, Args... args) {
            auto io_obj_executor = obj.get_executor();
            executor_info exec_info{};
            auto res = create_initial_netresult<R>();
            invoke_polyfill(
                fn,
                obj,
                std::forward<Args>(args)...,
                *res.diag,
                as_network_result<R>(res, create_tracker_executor(obj.get_executor(), &exec_info))
            );
            run_until_completion(io_obj_executor);
            BOOST_TEST(get_executor_info(io_obj_executor).num_posts > 0u);
            BOOST_TEST(exec_info.num_dispatches > 0u);
            return res;
        };
    }

    template <class Pfn>
    static signature async_noerrinfo(Pfn fn)
    {
        return [fn](IOObject& obj, Args... args) {
            auto io_obj_executor = obj.get_executor();
            executor_info exec_info{};
            auto res = create_initial_netresult<R>(false);
            invoke_polyfill(
                fn,
                obj,
                std::forward<Args>(args)...,
                as_network_result<R>(res, create_tracker_executor(obj.get_executor(), &exec_info))
            );
            run_until_completion(io_obj_executor);
            BOOST_TEST(get_executor_info(io_obj_executor).num_posts > 0u);
            BOOST_TEST(exec_info.num_dispatches > 0u);
            return res;
        };
    }

    // Used by channel functions
    template <class Pfn>
    static signature sync_errc_noerrinfo(Pfn fn)
    {
        return [fn](IOObject& obj, Args... args) {
            auto res = create_initial_netresult<R>(false);
            invoke_and_assign(res, fn, obj, std::forward<Args>(args)..., res.err);
            return res;
        };
    }
};

template <class R, class Obj, class... Args>
class netfun_maker_mem
{
    using impl = netfun_maker_impl<R, Obj&, Args...>;

public:
    using signature = std::function<network_result<R>(Obj&, Args...)>;
    using sig_sync_errc = R (Obj::*)(Args..., error_code&, diagnostics&);
    using sig_sync_errc_noerrinfo = R (Obj::*)(Args..., error_code&);
    using sig_sync_exc = R (Obj::*)(Args...);
    using sig_async_errinfo = void (Obj::*)(Args..., diagnostics&, as_network_result<R>&&);
    using sig_async_noerrinfo = void (Obj::*)(Args..., as_network_result<R>&&);

    static signature sync_errc(sig_sync_errc pfn) { return impl::sync_errc(pfn); }
    static signature sync_errc_noerrinfo(sig_sync_errc_noerrinfo pfn)
    {
        return impl::sync_errc_noerrinfo(pfn);
    }
    static signature sync_exc(sig_sync_exc pfn) { return impl::sync_exc(pfn); }
    static signature async_errinfo(sig_async_errinfo pfn) { return impl::async_errinfo(pfn); }
    static signature async_noerrinfo(sig_async_noerrinfo pfn) { return impl::async_noerrinfo(pfn); }
};

template <class R, class... Args>
class netfun_maker_fn
{
    using impl = netfun_maker_impl<R, Args...>;

public:
    using signature = std::function<network_result<R>(Args...)>;
    using sig_sync_errc = R (*)(Args..., error_code&, diagnostics&);
    using sig_sync_errc_noerrinfo = R (*)(Args..., error_code&);
    using sig_sync_exc = R (*)(Args...);
    using sig_async_errinfo = void (*)(Args..., diagnostics&, as_network_result<R>&&);
    using sig_async_noerrinfo = void (*)(Args..., as_network_result<R>&&);

    static signature sync_errc(sig_sync_errc pfn) { return impl::sync_errc(pfn); }
    static signature sync_errc_noerrinfo(sig_sync_errc_noerrinfo pfn)
    {
        return impl::sync_errc_noerrinfo(pfn);
    }
    static signature sync_exc(sig_sync_exc pfn) { return impl::sync_exc(pfn); }
    static signature async_errinfo(sig_async_errinfo pfn) { return impl::async_errinfo(pfn); }
    static signature async_noerrinfo(sig_async_noerrinfo pfn) { return impl::async_noerrinfo(pfn); }
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
