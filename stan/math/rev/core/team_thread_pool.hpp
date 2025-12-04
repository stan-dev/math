#ifndef STAN_MATH_REV_CORE_TEAM_THREAD_POOL_HPP
#define STAN_MATH_REV_CORE_TEAM_THREAD_POOL_HPP

#include <stan/math/rev/core/chainablestack.hpp>

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace stan {
namespace math {

/**
 * Team (epoch) thread pool for low-overhead parallel regions.
 *
 * - Creates (hw-1) worker threads once.
 * - Caller participates with tid=0.
 * - parallel_region(n, fn): runs fn(tid) for tid in [0, n).
 * - Nested parallelism: if called from a worker thread, runs serial.
 *
 * Designed for reduce_sum/map_rect style internal parallelism.
 */
class TeamThreadPool {
 public:
  static TeamThreadPool& instance() {
    static TeamThreadPool pool;
    return pool;
  }

  TeamThreadPool(const TeamThreadPool&) = delete;
  TeamThreadPool& operator=(const TeamThreadPool&) = delete;

  // Number of worker threads (excluding caller)
  std::size_t worker_count() const noexcept { return workers_.size(); }

  // Total participants available = worker_count + 1 (caller)
  std::size_t team_size() const noexcept { return workers_.size() + 1; }

  template <typename F>
  void parallel_region(std::size_t n, F&& fn) {
    if (n == 0) return;

    // If called from a worker, run serial to avoid nested deadlocks.
    if (in_worker_) {
      fn(std::size_t{0});
      return;
    }

    const std::size_t max_team = team_size();
    if (max_team == 1) {
      fn(std::size_t{0});
      return;
    }
    if (n > max_team) n = max_team;
    if (n == 1) {
      fn(std::size_t{0});
      return;
    }

    // Stable storage for callable during this region
    using Fn = std::decay_t<F>;
    Fn fn_copy = std::forward<F>(fn);

    // Publish region
    remaining_.store(n - 1, std::memory_order_release);  // workers only
    region_n_.store(n, std::memory_order_release);
    region_ctx_.store(static_cast<void*>(&fn_copy), std::memory_order_release);
    region_call_.store(&call_impl<Fn>, std::memory_order_release);

    epoch_.fetch_add(1, std::memory_order_acq_rel);

    // Wake workers
    {
      std::lock_guard<std::mutex> lk(wake_m_);
    }
    wake_cv_.notify_all();

    // Caller participates as tid=0
    in_worker_ = true;
    fn_copy(0);
    in_worker_ = false;

    // Wait for workers 1..n-1
    std::unique_lock<std::mutex> lk(done_m_);
    done_cv_.wait(lk, [&] {
      return remaining_.load(std::memory_order_acquire) == 0;
    });
  }

 private:
  using call_fn_t = void (*)(void*, std::size_t);

  template <typename Fn>
  static void call_impl(void* ctx, std::size_t tid) {
    (*static_cast<Fn*>(ctx))(tid);
  }

  TeamThreadPool()
      : stop_(false), epoch_(0), region_n_(0), region_ctx_(nullptr),
        region_call_(nullptr), remaining_(0) {
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) hw = 2;

    // hw-1 worker threads; caller is +1 participant.
    const unsigned num_workers = (hw > 1) ? (hw - 1) : 1;

    workers_.reserve(num_workers);
    for (unsigned i = 0; i < num_workers; ++i) {
      const std::size_t tid = static_cast<std::size_t>(i + 1);  // workers are 1..N
      workers_.emplace_back([this, tid] {
        // Per-worker AD tape initialized once
        static thread_local ChainableStack ad_tape;
        in_worker_ = true;

        std::size_t seen = epoch_.load(std::memory_order_acquire);
        for (;;) {
          // Sleep until epoch changes or stop requested
          {
            std::unique_lock<std::mutex> lk(wake_m_);
            wake_cv_.wait(lk, [&] {
              return stop_.load(std::memory_order_acquire)
                     || epoch_.load(std::memory_order_acquire) != seen;
            });
          }
          if (stop_.load(std::memory_order_acquire)) break;

          const std::size_t e = epoch_.load(std::memory_order_acquire);
          seen = e;

          const std::size_t n = region_n_.load(std::memory_order_acquire);
          if (tid >= n) {
            continue;  // not participating this region
          }

          void* ctx = region_ctx_.load(std::memory_order_acquire);
          call_fn_t call = region_call_.load(std::memory_order_acquire);
          if (call) {
            call(ctx, tid);
          }

          if (remaining_.fetch_sub(1, std::memory_order_acq_rel) == 1) {
            std::lock_guard<std::mutex> lk(done_m_);
            done_cv_.notify_one();
          }
        }

        in_worker_ = false;
      });
    }
  }

  ~TeamThreadPool() {
    stop_.store(true, std::memory_order_release);
    {
      std::lock_guard<std::mutex> lk(wake_m_);
    }
    wake_cv_.notify_all();
    for (auto& t : workers_) {
      if (t.joinable()) t.join();
    }
  }

  static inline thread_local bool in_worker_ = false;

  std::vector<std::thread> workers_;
  std::atomic<bool> stop_;

  // Region publication
  std::atomic<std::size_t> epoch_;
  std::atomic<std::size_t> region_n_;
  std::atomic<void*> region_ctx_;
  std::atomic<call_fn_t> region_call_;

  // Worker wake
  std::mutex wake_m_;
  std::condition_variable wake_cv_;

  // Completion
  std::atomic<std::size_t> remaining_;
  std::mutex done_m_;
  std::condition_variable done_cv_;
};

}  // namespace math
}  // namespace stan

#endif
