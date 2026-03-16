#pragma once

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <exception>
#include <functional>
#include <iostream>
#include <mutex>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

// Stan include (adjust include path if needed)
#include <stan/math/rev/core/chainablestack.hpp>

namespace stan {
namespace math {

/**
 * TeamThreadPool
 *
 * A "team" pool runs *collective parallel regions*:
 *  - parallel_region(n, fn) runs fn(tid) for tid=0..n-1
 *  - tid=0 runs on the caller thread
 *  - tid=1..n-1 run on persistent worker threads
 *
 * This is optimized for repeated short parallel regions (like
 * reduce_sum/map_rect), avoiding task-queue overhead.
 *
 * This version uses a single atomic wake generation counter (wake_gen_) and
 * removes the older "epoch" concept entirely.
 */
class TeamThreadPool {
 public:
  TeamThreadPool(const TeamThreadPool&) = delete;
  TeamThreadPool& operator=(const TeamThreadPool&) = delete;

  static TeamThreadPool& instance() {
    static TeamThreadPool pool;
    return pool;
  }

  /** Set total participants INCLUDING caller (tid=0). Call before instance().
   */
  static void set_num_threads(int n) {
    if (n < 1)
      n = 1;
    configured_threads_.store(n, std::memory_order_release);
  }

  /** Returns configured total participants (caller + workers). */
  static int get_num_threads() {
    return configured_threads_.load(std::memory_order_acquire);
  }

  /** Total participants INCLUDING caller (tid=0). */
  std::size_t team_size() const noexcept { return workers_.size() + 1; }

  /** Number of worker threads (excludes caller). */
  std::size_t worker_count() const noexcept { return workers_.size(); }

  template <typename F>
  void parallel_region(std::size_t n, F&& fn) {
    if (n == 0)
      return;

    // Clamp to actual team size
    const std::size_t max_team = team_size();
    if (max_team <= 1) {
      fn(std::size_t{0});
      return;
    }
    if (n > max_team)
      n = max_team;
    if (n <= 1) {
      fn(std::size_t{0});
      return;
    }

    // Nested parallelism: execute ALL tids serially on this worker thread.
    // This preserves correctness (every tid's chunk runs).
    if (in_worker_) {
      std::exception_ptr ep = nullptr;
      for (std::size_t tid = 0; tid < n; ++tid) {
        try {
          fn(tid);
        } catch (...) {
          if (!ep)
            ep = std::current_exception();
        }
      }
      if (ep)
        std::rethrow_exception(ep);
      return;
    }

    // Only one active region at a time.
    std::unique_lock<std::mutex> region_lock(region_m_);

    using Fn = std::decay_t<F>;
    Fn fn_copy = std::forward<F>(fn);

    // Exception propagation (first exception wins)
    std::exception_ptr eptr = nullptr;
    {
      std::lock_guard<std::mutex> lk(exc_m_);
      exc_ptr_ = &eptr;
    }

    // Publish region state (workers read with acquire)
    remaining_.store(n - 1, std::memory_order_release);  // workers only
    region_n_.store(n, std::memory_order_release);
    region_ctx_.store(static_cast<void*>(&fn_copy), std::memory_order_release);
    region_call_.store(&call_impl<Fn>, std::memory_order_release);

    // Start region: bump wake generation under wake_m_
    {
      std::lock_guard<std::mutex> lk(wake_m_);
      wake_gen_.fetch_add(1, std::memory_order_release);
    }
    wake_cv_.notify_all();

    // Caller participates as tid=0
    in_worker_ = true;
    try {
      fn_copy(0);
    } catch (...) {
      std::lock_guard<std::mutex> lk(exc_m_);
      if (!eptr)
        eptr = std::current_exception();
    }
    in_worker_ = false;

    // Wait for workers 1..n-1
    {
      std::unique_lock<std::mutex> lk(done_m_);
      done_cv_.wait(
          lk, [&] { return remaining_.load(std::memory_order_acquire) == 0; });
    }

    // Hygiene: deactivate region state
    region_call_.store(nullptr, std::memory_order_release);
    region_ctx_.store(nullptr, std::memory_order_release);
    region_n_.store(0, std::memory_order_release);

    // Clear exception slot
    {
      std::lock_guard<std::mutex> lk(exc_m_);
      exc_ptr_ = nullptr;
    }

    if (eptr)
      std::rethrow_exception(eptr);
  }

  static bool in_worker_thread() noexcept { return in_worker_; }

 private:
  using call_fn_t = void (*)(void*, std::size_t);

  template <typename Fn>
  static void call_impl(void* ctx, std::size_t tid) {
    (*static_cast<Fn*>(ctx))(tid);
  }

  static std::size_t configured_cap_(std::size_t hw) {
    int cfg = configured_threads_.load(std::memory_order_acquire);
    std::size_t cap = (cfg > 0) ? static_cast<std::size_t>(cfg) : hw;
    if (cap < 1)
      cap = 1;
    if (cap > hw)
      cap = hw;  // don't exceed hardware threads by default
    return cap;
  }

  TeamThreadPool()
      : stop_(false),
        wake_gen_(0),
        region_n_(0),
        region_ctx_(nullptr),
        region_call_(nullptr),
        remaining_(0),
        exc_ptr_(nullptr),
        ready_count_(0) {
    unsigned hw_u = std::thread::hardware_concurrency();
    if (hw_u == 0)
      hw_u = 2;
    const std::size_t hw = static_cast<std::size_t>(hw_u);

    const std::size_t cap = configured_cap_(hw);
    const std::size_t num_workers = (cap > 1) ? (cap - 1) : 0;

    std::cerr << "[TeamThreadPool(ctor)] cap=" << cap
              << " (workers=" << num_workers << ") hw=" << hw << "\n";

    workers_.reserve(num_workers);
    for (std::size_t i = 0; i < num_workers; ++i) {
      const std::size_t tid = i + 1;  // workers are 1..num_workers

      workers_.emplace_back([this, tid] {
        // Per-worker AD tape initialized once.
        static thread_local ChainableStack ad_tape;
        (void)ad_tape;

        in_worker_ = true;

        // Startup barrier: ensure each worker reached the wait loop at least
        // once.
        {
          std::lock_guard<std::mutex> lk(wake_m_);
          ready_count_.fetch_add(1, std::memory_order_release);
        }
        ready_cv_.notify_one();

        // Track wake generation.
        std::size_t seen_gen = wake_gen_.load(std::memory_order_acquire);

        while (!stop_.load(std::memory_order_acquire)) {
          // Wait until wake_gen changes (or stop)
          {
            std::unique_lock<std::mutex> lk(wake_m_);
            wake_cv_.wait(lk, [&] {
              return stop_.load(std::memory_order_acquire)
                     || wake_gen_.load(std::memory_order_acquire) != seen_gen;
            });
            if (stop_.load(std::memory_order_acquire))
              break;

            // IMPORTANT: update while holding wake_m_ so we can't miss rapid
            // increments
            seen_gen = wake_gen_.load(std::memory_order_acquire);
          }

          const std::size_t n = region_n_.load(std::memory_order_acquire);
          if (tid >= n) {
            continue;  // not participating in this region
          }

          // Always decrement once for participating workers.
          struct DoneGuard {
            std::atomic<std::size_t>& rem;
            std::mutex& m;
            std::condition_variable& cv;
            ~DoneGuard() {
              if (rem.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                std::lock_guard<std::mutex> lk(m);
                cv.notify_one();
              }
            }
          } guard{remaining_, done_m_, done_cv_};

          void* ctx = region_ctx_.load(std::memory_order_acquire);
          call_fn_t call = region_call_.load(std::memory_order_acquire);

          if (call) {
            try {
              call(ctx, tid);
            } catch (...) {
              std::lock_guard<std::mutex> lk(exc_m_);
              if (exc_ptr_ && *exc_ptr_ == nullptr) {
                *exc_ptr_ = std::current_exception();
              }
            }
          }
        }

        in_worker_ = false;
      });
    }

    // Wait for workers to reach the wait loop once before returning.
    {
      std::unique_lock<std::mutex> lk(wake_m_);
      ready_cv_.wait(lk, [&] {
        return ready_count_.load(std::memory_order_acquire) == workers_.size();
      });
    }

    std::cerr << "[TeamThreadPool(ctor)] ready workers=" << workers_.size()
              << " team_size=" << team_size() << "\n";
  }

  ~TeamThreadPool() {
    stop_.store(true, std::memory_order_release);
    {
      std::lock_guard<std::mutex> lk(wake_m_);
      wake_gen_.fetch_add(1, std::memory_order_release);
    }
    wake_cv_.notify_all();
    for (auto& th : workers_) {
      if (th.joinable())
        th.join();
    }
  }

  // ---- configuration ----
  static inline std::atomic<int> configured_threads_{0};

  // ---- worker threads ----
  std::vector<std::thread> workers_;
  static inline thread_local bool in_worker_ = false;

  // ---- region serialization ----
  std::mutex region_m_;

  // ---- wake mechanism (generation counter) ----
  std::mutex wake_m_;
  std::condition_variable wake_cv_;
  std::atomic<std::size_t> wake_gen_;
  std::atomic<bool> stop_;

  // ---- region state (published by caller, read by workers) ----
  std::atomic<std::size_t> region_n_;
  std::atomic<void*> region_ctx_;
  std::atomic<call_fn_t> region_call_;

  // ---- completion ----
  std::atomic<std::size_t> remaining_;
  std::mutex done_m_;
  std::condition_variable done_cv_;

  // ---- exception propagation ----
  std::mutex exc_m_;
  std::exception_ptr* exc_ptr_;

  // ---- startup barrier ----
  std::atomic<std::size_t> ready_count_;
  std::condition_variable ready_cv_;
};

}  // namespace math
}  // namespace stan
