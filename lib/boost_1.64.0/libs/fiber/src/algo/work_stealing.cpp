
//          Copyright Oliver Kowalke 2015.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
//

#include "boost/fiber/algo/work_stealing.hpp"

#include <random>

#include <boost/assert.hpp>

#include "boost/fiber/type.hpp"

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

namespace boost {
namespace fibers {
namespace algo {

void
work_stealing::init_( std::size_t max_idx) {
    schedulers_.resize( max_idx + 1);
}

work_stealing::work_stealing( std::size_t max_idx, std::size_t idx, bool suspend) :
    idx_{ idx },
    max_idx_{ max_idx },
    suspend_{ suspend } {
    static std::once_flag flag;
    std::call_once( flag, & work_stealing::init_, max_idx_);
    schedulers_[idx_] = this;

}

void
work_stealing::awakened( context * ctx) noexcept {
    if ( ! ctx->is_context( type::pinned_context) ) {
        ctx->detach();
    }
    rqueue_.push( ctx);
}

context *
work_stealing::pick_next() noexcept {
    context * ctx = rqueue_.pop();
    if ( nullptr != ctx) {
        if ( ! ctx->is_context( type::pinned_context) ) {
            context::active()->attach( ctx);
        }
    } else {
        static thread_local std::minstd_rand generator;
        static std::uniform_int_distribution< std::size_t > distribution{ 0, max_idx_ };
        std::size_t idx = 0;
        do {
            idx = distribution( generator);
        } while ( idx == idx_);
        ctx = schedulers_[idx]->steal();
        if ( nullptr != ctx) {
            BOOST_ASSERT( ! ctx->is_context( type::pinned_context) );
            context::active()->attach( ctx);
        }
    }
    return ctx;
}

void
work_stealing::suspend_until( std::chrono::steady_clock::time_point const& time_point) noexcept {
    if ( suspend_) {
        if ( (std::chrono::steady_clock::time_point::max)() == time_point) {
            std::unique_lock< std::mutex > lk{ mtx_ };
            cnd_.wait( lk, [this](){ return flag_; });
            flag_ = false;
        } else {
            std::unique_lock< std::mutex > lk{ mtx_ };
            cnd_.wait_until( lk, time_point, [this](){ return flag_; });
            flag_ = false;
        }
    }
}

void
work_stealing::notify() noexcept {
    if ( suspend_) {
        std::unique_lock< std::mutex > lk{ mtx_ };
        flag_ = true;
        lk.unlock();
        cnd_.notify_all();
    }
}

std::vector< work_stealing * > work_stealing::schedulers_{};

}}}

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif
