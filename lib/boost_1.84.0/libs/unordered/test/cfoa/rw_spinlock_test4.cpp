// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/unordered/detail/foa/rw_spinlock.hpp>
#include <boost/compat/shared_lock.hpp>
#include <mutex>

using boost::unordered::detail::foa::rw_spinlock;

static rw_spinlock sp;
static rw_spinlock sp2;

int main()
{
    sp.lock();
    sp2.lock_shared();
    sp2.lock_shared();

    sp.unlock();
    sp2.unlock_shared();
    sp2.unlock_shared();

    {
        std::lock_guard<rw_spinlock> lock( sp );
        boost::compat::shared_lock<rw_spinlock> lock2( sp2 );
        boost::compat::shared_lock<rw_spinlock> lock3( sp2 );
    }
}
