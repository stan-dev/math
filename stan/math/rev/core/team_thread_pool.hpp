#ifndef STAN_MATH_REV_CORE_TEAM_THREAD_POOL_HPP
#define STAN_MATH_REV_CORE_TEAM_THREAD_POOL_HPP

#include <stan/math/rev/core/chainablestack.hpp>

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <mutex>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace stan {
namespace math {

/**
 * Team (epoch) thread pool for low-overhead parallel regions.
 *
 * - Workers are created once.
 * - Caller participates as tid=0.
 * - parallel_region(n, fn): runs fn(tid) for tid in [0, n).
 * - Nested parallelism: if called from a worker thread, runs serial.
 * - set_num_threads(k) must be called before instance() to size the pool.
 */
class TeamThreadPool {
 public:
  // Call before first instance() to control pool size.
  // Meaning: total participants INCLUDING caller (tid=0).
  static void set_num_threads(std::size_t n) noexcept {
    if (n < 1) n = 1;
    user_cap_().store(n, std::memory_order_release);
  }

  static std::size_t get_num_threads() noexcept {
    return user_cap_().load(std::memory_order_acquire);
  }

  static TeamThreadPool& instance() {
    static TeamThreadPool pool;
    return pool;
  }

  // Worker threads (excluding caller)
  std::size_t worker_count() const noexcept { return workers_.size(); }

  // Total possible participants INCLUDING caller
  std::size_t team_size() const noexcept { return workers_.size() + 1; }

  template <typename F>
  void parallel_region(std::size_t n, F&& fn) {
    if (n == 0) return;

    // Nested parallelism guard: if already on a worker, run serial.
    if (in_worker_) {
      fn(std::size_t{0});
      return;
    }

    // Only one active region at a time (this is required for a single shared epoch design).
    std::unique_lock<std::mutex> region_lock(region_m_);

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

    using Fn = std::decay_t<F>;
    Fn fn_copy = std::forward<F>(fn);  // stable storage during this call

    // Exception propagation: capture first exception from any participant.
    std::exception_ptr eptr = nullptr;
    {
      std::lock_guard<std::mutex> lk(exc_m_);
      exc_ptr_ = &eptr;
    }

    // Publish region state BEFORE bumping epoch.
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
    try {
      fn_copy(0);
    } catch (...) {
      std::lock_guard<std::mutex> lk(exc_m_);
      if (eptr == nullptr) eptr = std::current_exception();
    }
    in_worker_ = false;

    // Wait for workers 1..n-1
    std::unique_lock<std::mutex> lk(done_m_);
    done_cv_.wait(lk, [&] {
      return remaining_.load(std::memory_order_acquire) == 0;
    });

    // Clear region participation; not strictly necessary but good hygiene.
    region_n_.store(0, std::memory_order_release);

    // Rethrow exception (if any)
    if (eptr) std::rethrow_exception(eptr);
  }

 private:
  using call_fn_t = void (*)(void*, std::size_t);

  template <typename Fn>
  static void call_impl(void* ctx, std::size_t tid) {
    (*static_cast<Fn*>(ctx))(tid);
  }

  // Function-local static avoids static initialization order issues.
  static std::atomic<std::size_t>& user_cap_() {
    static std::atomic<std::size_t> cap{0};  // 0 => unset
    return cap;
  }

  static std::size_t env_num_threads_() noexcept {
    const char* s = std::getenv("STAN_NUM_THREADS");
    if (!s || !*s) return 0;
    char* end = nullptr;
    long v = std::strtol(s, &end, 10);
    if (end == s || v <= 0) return 0;
    return static_cast<std::size_t>(v);
  }

  static std::size_t configured_cap_(std::size_t hw) noexcept {
    // Priority: explicit set_num_threads > STAN_NUM_THREADS > hw
    std::size_t cap = user_cap_().load(std::memory_order_acquire);
    if (cap == 0) cap = env_num_threads_();
    if (cap == 0) cap = hw;

    if (cap < 1) cap = 1;
    if (cap > hw) cap = hw;  // don’t oversubscribe by default
    return cap;
  }

  TeamThreadPool()
      : stop_(false),
        epoch_(0),
        region_n_(0),
        region_ctx_(nullptr),
        region_call_(nullptr),
        remaining_(0),
        exc_ptr_(nullptr) {
    unsigned hw_u = std::thread::hardware_concurrency();
    if (hw_u == 0) hw_u = 2;
    const std::size_t hw = static_cast<std::size_t>(hw_u);

    // Total participants includes caller.
    const std::size_t cap = configured_cap_(hw);
    const std::size_t num_workers = (cap > 1) ? (cap - 1) : 0;

    workers_.reserve(num_workers);
    for (std::size_t i = 0; i < num_workers; ++i) {
      const std::size_t tid = i + 1;  // workers are 1..num_workers
      workers_.emplace_back([this, tid] {
        static thread_local ChainableStack ad_tape;
        in_worker_ = true;

        std::size_t seen = epoch_.load(std::memory_order_acquire);
        for (;;) {
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
          if (tid >= n) continue;  // not participating this region

          // Ensure we ALWAYS decrement remaining_ once for participating workers.
          struct DoneGuard {
            std::atomic<std::size_t>& rem;
            std::mutex& m;
            std::condition_variable& cv;
            bool active{true};
            ~DoneGuard() {
              if (!active) return;
              if (rem.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                std::lock_guard<std::mutex> lk(m);
                cv.notify_one();
              }
            }
          } guard{remaining_, done_m_, done_cv_};

          void* ctx = region_ctx_.load(std::memory_order_acquire);
          call_fn_t call = region_call_.load(std::memory_order_acquire);
          if (!call) continue;

          try {
            call(ctx, tid);
          } catch (...) {
            std::lock_guard<std::mutex> lk(exc_m_);
            if (exc_ptr_ && *exc_ptr_ == nullptr) {
              *exc_ptr_ = std::current_exception();
            }
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

  // Serialize regions (single shared-region design)
  std::mutex region_m_;

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

  // Exception plumbing
  std::mutex exc_m_;
  std::exception_ptr* exc_ptr_;
};

}  // namespace math
}  // namespace stan

#endif
