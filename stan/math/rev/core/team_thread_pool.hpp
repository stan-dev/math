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
  // Call this before first use of TeamThreadPool::instance()
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

  TeamThreadPool(const TeamThreadPool&) = delete;
  TeamThreadPool& operator=(const TeamThreadPool&) = delete;

  // Number of worker threads (excluding caller)
  std::size_t worker_count() const noexcept { return workers_.size(); }

  // Total participants available = worker_count + 1 (caller)
  std::size_t team_size() const noexcept { return workers_.size() + 1; }

  template <typename F>
  void parallel_region(std::size_t n, F&& fn) {
    //std::cout << "#################### parallel_region, n = " << n << std::endl;
    if (n == 0) return;

    // If called from a worker, run serial to avoid nested deadlocks.
    //std::cout << "in_worker_ = " << in_worker_ << std::endl;
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

    //std::cout << "waiting for workers" << std::endl;
    // Wait for workers 1..n-1
    std::unique_lock<std::mutex> lk(done_m_);
    done_cv_.wait(lk, [&] {
      return remaining_.load(std::memory_order_acquire) == 0;
    });
    //std::cout << "#################### done" << std::endl << std::endl;
  }

private:
  // Function-local static avoids static init order fiasco.
  static std::atomic<std::size_t>& user_cap_() {
    static std::atomic<std::size_t> cap{0};  // 0 => "unset"
    return cap;
  }

  static std::size_t configured_cap_(std::size_t hw) {
    // priority: user cap > env var > hw
    std::size_t cap = user_cap_().load(std::memory_order_acquire);
    if (cap == 0) {
      cap = env_num_threads_();   // if you have STAN_NUM_THREADS support
    }
    if (cap == 0) cap = hw;

    if (cap < 1) cap = 1;
    if (cap > hw) cap = hw;       // prevent oversubscription by default
    return cap;
  }

  
  using call_fn_t = void (*)(void*, std::size_t);

  template <typename Fn>
  static void call_impl(void* ctx, std::size_t tid) {
    (*static_cast<Fn*>(ctx))(tid);
  }

  static size_t env_num_threads_() {
    size_t num_threads = 1;
#ifdef STAN_THREADS
    const char* env_stan_num_threads = std::getenv("STAN_NUM_THREADS");
    if (env_stan_num_threads != nullptr) {
      try {
	const int env_num_threads
          = boost::lexical_cast<int>(env_stan_num_threads);
	if (env_num_threads > 0) {
	  num_threads = env_num_threads;
	} else if (env_num_threads == -1) {
	  num_threads = std::thread::hardware_concurrency();
	} else {
	  invalid_argument("get_num_threads(int)", "STAN_NUM_THREADS",
			   env_stan_num_threads,
			   "The STAN_NUM_THREADS environment variable is '",
			   "' but it must be positive or -1");
	}
      } catch (const boost::bad_lexical_cast&) {
	invalid_argument("get_num_threads(int)", "STAN_NUM_THREADS",
			 env_stan_num_threads,
			 "The STAN_NUM_THREADS environment variable is '",
			 "' but it must be a positive number or -1");
      }
    }
#endif
  return num_threads;
}


  TeamThreadPool()
      : stop_(false), epoch_(0), region_n_(0), region_ctx_(nullptr),
        region_call_(nullptr), remaining_(0) {

    unsigned hw_u = std::thread::hardware_concurrency();
    if (hw_u == 0) hw_u = 2;
    const std::size_t hw = static_cast<std::size_t>(hw_u);
    
    const std::size_t cap = configured_cap_(hw);
    const std::size_t num_workers = (cap > 1) ? (cap - 1) : 0;

    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl
	      << "hw = " << hw << std::endl 
	      << "num_workers = " << num_workers << std::endl
	      << "cap = " << cap << std::endl
	      << std::endl << std::endl;
    
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
    std::cout << "done with constructor" << std::endl;
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
