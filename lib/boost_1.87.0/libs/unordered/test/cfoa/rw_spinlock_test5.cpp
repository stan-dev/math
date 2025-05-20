// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/unordered/detail/foa/rw_spinlock.hpp>
#include <boost/compat/shared_lock.hpp>
#include <boost/core/lightweight_test.hpp>
#include <mutex>

using boost::unordered::detail::foa::rw_spinlock;

static rw_spinlock sp;

int main()
{
    {
        BOOST_TEST( sp.try_lock_shared() );
        BOOST_TEST( sp.try_lock_shared() );
        sp.unlock_shared();
        sp.unlock_shared();
    }

    {
        BOOST_TEST( sp.try_lock() );
        BOOST_TEST( !sp.try_lock_shared() );
        sp.unlock();
    }

    {
        std::lock_guard<rw_spinlock> lock( sp );
        BOOST_TEST( !sp.try_lock_shared() );
    }

    {
        boost::compat::shared_lock<rw_spinlock> lock( sp );
        BOOST_TEST( !sp.try_lock() );
        BOOST_TEST( sp.try_lock_shared() );
        sp.unlock_shared();
    }

    return boost::report_errors();
}
