// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/unordered/detail/foa/rw_spinlock.hpp>
#include <boost/core/lightweight_test.hpp>
#include <mutex>

using boost::unordered::detail::foa::rw_spinlock;

static rw_spinlock sp;
static rw_spinlock sp2;

int main()
{
    BOOST_TEST( sp.try_lock() );
    BOOST_TEST( !sp.try_lock() );
    BOOST_TEST( sp2.try_lock() );
    BOOST_TEST( !sp.try_lock() );
    BOOST_TEST( !sp2.try_lock() );
    sp.unlock();
    sp2.unlock();

    sp.lock();
    BOOST_TEST( !sp.try_lock() );
    sp2.lock();
    BOOST_TEST( !sp.try_lock() );
    BOOST_TEST( !sp2.try_lock() );
    sp.unlock();
    sp2.unlock();

    {
        std::lock_guard<rw_spinlock> lock( sp );
        BOOST_TEST( !sp.try_lock() );
        std::lock_guard<rw_spinlock> lock2( sp2 );
        BOOST_TEST( !sp.try_lock() );
        BOOST_TEST( !sp2.try_lock() );
    }

    return boost::report_errors();
}
