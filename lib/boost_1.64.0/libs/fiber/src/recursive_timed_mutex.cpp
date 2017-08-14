
//          Copyright Oliver Kowalke 2013.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "boost/fiber/recursive_timed_mutex.hpp"

#include <algorithm>
#include <functional>

#include "boost/fiber/exceptions.hpp"
#include "boost/fiber/scheduler.hpp"

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

namespace boost {
namespace fibers {

bool
recursive_timed_mutex::try_lock_until_( std::chrono::steady_clock::time_point const& timeout_time) noexcept {
    if ( std::chrono::steady_clock::now() > timeout_time) {
        return false;
    }
    context * active_ctx = context::active();
    // store this fiber in order to be notified later
    detail::spinlock_lock lk{ wait_queue_splk_ };
    if ( active_ctx == owner_) {
        ++count_;
        return true;
    } else if ( nullptr == owner_) {
        owner_ = active_ctx;
        count_ = 1;
        return true;
    }
    BOOST_ASSERT( ! active_ctx->wait_is_linked() );
    active_ctx->wait_link( wait_queue_);
    // suspend this fiber until notified or timed-out
    if ( ! active_ctx->wait_until( timeout_time, lk) ) {
        // remove fiber from wait-queue 
        lk.lock();
        wait_queue_.remove( * active_ctx);
        return false;
    }
    BOOST_ASSERT( ! active_ctx->wait_is_linked() );
    return active_ctx == owner_;
}

void
recursive_timed_mutex::lock() {
    context * active_ctx = context::active();
    // store this fiber in order to be notified later
    detail::spinlock_lock lk{ wait_queue_splk_ };
    if ( active_ctx == owner_) {
        ++count_;
        return;
    } else if ( nullptr == owner_) {
        owner_ = active_ctx;
        count_ = 1;
        return;
    }
    BOOST_ASSERT( ! active_ctx->wait_is_linked() );
    active_ctx->wait_link( wait_queue_);
    // suspend this fiber
    active_ctx->suspend( lk);
    BOOST_ASSERT( ! active_ctx->wait_is_linked() );
}

bool
recursive_timed_mutex::try_lock() noexcept {
    context * active_ctx = context::active();
    detail::spinlock_lock lk{ wait_queue_splk_ };
    if ( nullptr == owner_) {
        owner_ = active_ctx;
        count_ = 1;
    } else if ( active_ctx == owner_) {
        ++count_;
    }
    lk.unlock();
    // let other fiber release the lock
    active_ctx->yield();
    return active_ctx == owner_;
}

void
recursive_timed_mutex::unlock() {
    context * active_ctx = context::active();
    detail::spinlock_lock lk{ wait_queue_splk_ };
    if ( active_ctx != owner_) {
        throw lock_error{
                std::make_error_code( std::errc::operation_not_permitted),
                "boost fiber: no  privilege to perform the operation" };
    }
    if ( 0 == --count_) {
        if ( ! wait_queue_.empty() ) {
            context * ctx = & wait_queue_.front();
            wait_queue_.pop_front();
            owner_ = ctx;
            count_ = 1;
            active_ctx->schedule( ctx);
        } else {
            owner_ = nullptr;
            return;
        }
    }
}

}}

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif
