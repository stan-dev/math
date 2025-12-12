#ifndef STAN_MATH_REV_CORE_TEAM_THREAD_POOL_HPP
#define STAN_MATH_REV_CORE_TEAM_THREAD_POOL_HPP

#include <stan/math/rev/core/chainablestack.hpp>

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdio>     // fprintf, fflush
#include <cstdlib>    // getenv, strtol
#include <exception>  // exception_ptr
#include <mutex>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

// Debug mode: periodic dumps while waiting.
// #define STAN_TEAM_POOL_DEBUG_WAIT 1

#if defined(STAN_TEAM_POOL_DEBUG_WAIT)
#include <chrono>
#endif

namespace stan {
namespace math {

/**
 * Team (epoch) thread pool for low-overhead parallel regions.
 *
 * Logical tids:
 * - caller thread: tid=0
 * - worker threads: tid=1..cap-1
 *
 * Key correctness points:
 * - Single shared "region" state => serialize parallel_region() with region_m_.
 * - Wake generation (wake_gen_) protected by wake_m_ to prevent missed wakeups.
 * - Startup barrier ensures all workers are "armed" on wake_gen_ before first use,
 *   preventing late-start workers from missing the first region.
 *
 * Debug fields per logical tid:
 * - alive[tid]    : 1 once the worker thread started (0 means never started / died early)
 * - seen[tid]     : last epoch observed by that worker (0 means never saw any region)
 * - exec[tid]     : epoch currently executing (0 means not executing region work)
 * - dec[tid]      : count of decrements performed by that worker
 */
class TeamThreadPool {
 public:
  // Total participants INCLUDING caller (tid=0). Call before instance().
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

  std::size_t worker_count() const noexcept { return workers_.size(); }
  std::size_t team_size() const noexcept { return workers_.size() + 1; }

  void dump_state(const char* tag = "TeamThreadPool") const {
    const auto epoch = epoch_.load(std::memory_order_acquire);
    const auto n = region_n_.load(std::memory_order_acquire);
    const auto rem = remaining_.load(std::memory_order_acquire);

    std::fprintf(stderr,
                 "\n[%s] epoch=%zu region_n=%zu remaining=%zu team_size=%zu wake_gen=%zu ready=%zu/%zu\n",
                 tag, epoch, n, rem, team_size(), wake_gen_snapshot_(),
                 ready_count_.load(std::memory_order_acquire),
                 workers_.size());
    std::fprintf(stderr, "  tid | alive |  seen |  exec | dec | note\n");
    std::fprintf(stderr, "  ----+-------+-------+-------+-----+-----------------------------\n");

    for (std::size_t tid = 1; tid < worker_state_size_; ++tid) {
      const unsigned alive = worker_alive_[tid].load(std::memory_order_acquire);
      const std::size_t seen = worker_seen_epoch_[tid].load(std::memory_order_acquire);
      const std::size_t exec = worker_exec_epoch_[tid].load(std::memory_order_acquire);
      const std::size_t dec = worker_decrement_count_[tid].load(std::memory_order_acquire);

      const char* note = "";
      if (!alive) {
        note = "NOT ALIVE";
      } else if (tid < n) {
        if (seen < epoch) note = "participating but hasn't seen epoch";
        else if (exec == epoch) note = "executing";
        else if (exec != 0) note = "executing (old epoch?)";
        else note = "idle/finished";
      } else {
        if (seen == epoch) note = "saw epoch but not participating";
        else note = "not participating";
      }

      std::fprintf(stderr, "  %3zu |   %u   | %5zu | %5zu | %3zu | %s\n",
                   tid, alive, seen, exec, dec, note);
    }
    std::fflush(stderr);
  }

  template <typename F>
  void parallel_region(std::size_t n, F&& fn) {
    if (n == 0) return;

    // Nested parallelism guard.
    if (in_worker_) {
      fn(std::size_t{0});
      return;
    }

    // Single shared region state => serialize launches.
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
    Fn fn_copy = std::forward<F>(fn);

    // Exception propagation (first exception wins).
    std::exception_ptr eptr = nullptr;
    {
      std::lock_guard<std::mutex> lk(exc_m_);
      exc_ptr_ = &eptr;
    }

    // Publish region state BEFORE bumping epoch.
    remaining_.store(n - 1, std::memory_order_release);
    region_n_.store(n, std::memory_order_release);
    region_ctx_.store(static_cast<void*>(&fn_copy), std::memory_order_release);
    region_call_.store(&call_impl<Fn>, std::memory_order_release);

    const std::size_t new_epoch =
        epoch_.fetch_add(1, std::memory_order_acq_rel) + 1;

    // std::fprintf(stderr,
    //              "\n[TeamThreadPool(launch)] epoch=%zu n=%zu expected_workers=%zu team_size=%zu\n",
    //              new_epoch, n, n - 1, team_size());
    // std::fflush(stderr);

    // Wake workers using wake generation (prevents missed wakeups).
    {
      std::lock_guard<std::mutex> lk(wake_m_);
      ++wake_gen_;
    }
    wake_cv_.notify_all();

    // Caller participates (tid=0).
    in_worker_ = true;
    try {
      fn_copy(0);
    } catch (...) {
      std::lock_guard<std::mutex> lk(exc_m_);
      if (eptr == nullptr) eptr = std::current_exception();
    }
    in_worker_ = false;

    // Wait for workers 1..n-1.
    std::unique_lock<std::mutex> lk(done_m_);
#if defined(STAN_TEAM_POOL_DEBUG_WAIT)
    auto last_dump = std::chrono::steady_clock::now();
    while (remaining_.load(std::memory_order_acquire) != 0) {
      done_cv_.wait_for(lk, std::chrono::milliseconds(250));
      const auto now = std::chrono::steady_clock::now();
      if (now - last_dump > std::chrono::seconds(2)
          && remaining_.load(std::memory_order_acquire) != 0) {
        std::fprintf(stderr,
                     "[TeamThreadPool] waiting too long for epoch=%zu (remaining=%zu)\n",
                     new_epoch, remaining_.load(std::memory_order_acquire));
        dump_state("TeamThreadPool(wait)");
        last_dump = now;
      }
    }
#else
    done_cv_.wait(lk, [&] {
      return remaining_.load(std::memory_order_acquire) == 0;
    });
#endif

    // Hygiene.
    region_n_.store(0, std::memory_order_release);

    if (eptr) std::rethrow_exception(eptr);
  }

 private:
  using call_fn_t = void (*)(void*, std::size_t);

  template <typename Fn>
  static void call_impl(void* ctx, std::size_t tid) {
    (*static_cast<Fn*>(ctx))(tid);
  }

  static std::atomic<std::size_t>& user_cap_() {
    static std::atomic<std::size_t> cap{0};
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
    std::size_t cap = user_cap_().load(std::memory_order_acquire);
    if (cap == 0) cap = env_num_threads_();
    if (cap == 0) cap = hw;
    if (cap < 1) cap = 1;
    if (cap > hw) cap = hw;
    return cap;
  }

  std::size_t wake_gen_snapshot_() const {
    std::lock_guard<std::mutex> lk(wake_m_);
    return wake_gen_;
  }

  TeamThreadPool()
      : stop_(false),
        epoch_(0),
        region_n_(0),
        region_ctx_(nullptr),
        region_call_(nullptr),
        remaining_(0),
        exc_ptr_(nullptr),
        wake_gen_(0),
        ready_count_(0) {
    unsigned hw_u = std::thread::hardware_concurrency();
    if (hw_u == 0) hw_u = 2;
    const std::size_t hw = static_cast<std::size_t>(hw_u);

    const std::size_t cap = configured_cap_(hw);
    const std::size_t num_workers = (cap > 1) ? (cap - 1) : 0;

    worker_state_size_ = cap;

    // raw arrays so atomics aren't moved
    worker_alive_.reset(new std::atomic<unsigned>[cap]);
    worker_seen_epoch_.reset(new std::atomic<std::size_t>[cap]);
    worker_exec_epoch_.reset(new std::atomic<std::size_t>[cap]);
    worker_decrement_count_.reset(new std::atomic<std::size_t>[cap]);

    for (std::size_t i = 0; i < cap; ++i) {
      worker_alive_[i].store(0u, std::memory_order_relaxed);
      worker_seen_epoch_[i].store(0, std::memory_order_relaxed);
      worker_exec_epoch_[i].store(0, std::memory_order_relaxed);
      worker_decrement_count_[i].store(0, std::memory_order_relaxed);
    }

    std::fprintf(stderr,
                 "[TeamThreadPool(ctor)] cap=%zu (workers=%zu) hw=%zu\n",
                 cap, num_workers, hw);
    std::fflush(stderr);

    workers_.reserve(num_workers);
    for (std::size_t i = 0; i < num_workers; ++i) {
      const std::size_t tid = i + 1;
      workers_.emplace_back([this, tid] {
        static thread_local ChainableStack ad_tape;
        in_worker_ = true;

        if (tid < worker_state_size_) {
          worker_alive_[tid].store(1u, std::memory_order_release);
        }

        // "Arm" this worker on the current wake generation so it can't miss
        // the first region wake.
        std::size_t local_gen;
        {
          std::unique_lock<std::mutex> lk(wake_m_);
          local_gen = wake_gen_;
        }

        // Signal readiness AFTER arming.
        ready_count_.fetch_add(1, std::memory_order_acq_rel);
        {
          std::lock_guard<std::mutex> lk(ready_m_);
        }
        ready_cv_.notify_one();

        for (;;) {
          // Wait for wake_gen_ to change (or stop_).
          {
            std::unique_lock<std::mutex> lk(wake_m_);
            wake_cv_.wait(lk, [&] {
              return stop_.load(std::memory_order_acquire)
                     || wake_gen_ != local_gen;
            });
            if (stop_.load(std::memory_order_acquire)) break;
            local_gen = wake_gen_;
          }

          // Observe epoch and region parameters after wake.
          const std::size_t e = epoch_.load(std::memory_order_acquire);

          if (tid < worker_state_size_) {
            worker_seen_epoch_[tid].store(e, std::memory_order_release);
          }

          const std::size_t n = region_n_.load(std::memory_order_acquire);
          if (tid >= n) {
            continue;
          }

          // Always decrement once for participating workers.
          struct DoneGuard {
            std::atomic<std::size_t>& rem;
            std::mutex& m;
            std::condition_variable& cv;
            std::atomic<std::size_t>& dec_count;
            ~DoneGuard() {
              dec_count.fetch_add(1, std::memory_order_relaxed);
              if (rem.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                std::lock_guard<std::mutex> lk(m);
                cv.notify_one();
              }
            }
          } guard{remaining_, done_m_, done_cv_, worker_decrement_count_[tid]};

          if (tid < worker_state_size_) {
            worker_exec_epoch_[tid].store(e, std::memory_order_release);
          }

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

          if (tid < worker_state_size_) {
            worker_exec_epoch_[tid].store(0, std::memory_order_release);
          }
        }

        in_worker_ = false;
      });
    }

    // Startup barrier: ensure all workers are armed and waiting-ready.
    {
      std::unique_lock<std::mutex> lk(ready_m_);
      ready_cv_.wait(lk, [&] {
        return ready_count_.load(std::memory_order_acquire) == workers_.size();
      });
    }
    std::fprintf(stderr, "[TeamThreadPool(ctor)] all workers ready: %zu\n",
                 workers_.size());
    std::fflush(stderr);
  }

  ~TeamThreadPool() {
    stop_.store(true, std::memory_order_release);
    {
      // bump wake_gen_ so workers wake and see stop_
      std::lock_guard<std::mutex> lk(wake_m_);
      ++wake_gen_;
    }
    wake_cv_.notify_all();

    for (auto& t : workers_) {
      if (t.joinable()) t.join();
    }
  }

  static inline thread_local bool in_worker_ = false;

  std::vector<std::thread> workers_;
  std::atomic<bool> stop_;

  // Serialize regions.
  std::mutex region_m_;

  // Region publication.
  std::atomic<std::size_t> epoch_;
  std::atomic<std::size_t> region_n_;
  std::atomic<void*> region_ctx_;
  std::atomic<call_fn_t> region_call_;

  // Wake workers (wake_gen_ is protected by wake_m_).
  mutable std::mutex wake_m_;
  std::condition_variable wake_cv_;
  std::size_t wake_gen_;

  // Startup barrier.
  std::mutex ready_m_;
  std::condition_variable ready_cv_;
  std::atomic<std::size_t> ready_count_;

  // Completion.
  std::atomic<std::size_t> remaining_;
  std::mutex done_m_;
  std::condition_variable done_cv_;

  // Exceptions.
  std::mutex exc_m_;
  std::exception_ptr* exc_ptr_;

  // Debug state (arrays of atomics to avoid std::vector<atomic<...>> move issues).
  std::size_t worker_state_size_{0};
  std::unique_ptr<std::atomic<unsigned>[]> worker_alive_;
  std::unique_ptr<std::atomic<std::size_t>[]> worker_seen_epoch_;
  std::unique_ptr<std::atomic<std::size_t>[]> worker_exec_epoch_;
  std::unique_ptr<std::atomic<std::size_t>[]> worker_decrement_count_;
};

}  // namespace math
}  // namespace stan

#endif
