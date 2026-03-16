#ifndef STAN_MATH_REV_CORE_SIMPLE_THREAD_POOL_HPP
#define STAN_MATH_REV_CORE_SIMPLE_THREAD_POOL_HPP

#include <stan/math/rev/core/chainablestack.hpp>

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace stan {
namespace math {

class SimpleThreadPool {
 public:
  static SimpleThreadPool& instance() {
    static SimpleThreadPool pool;
    return pool;
  }

  SimpleThreadPool(const SimpleThreadPool&) = delete;
  SimpleThreadPool& operator=(const SimpleThreadPool&) = delete;

  std::size_t thread_count() const noexcept { return workers_.size(); }

  template <typename F, typename... Args>
  auto submit(F&& f, Args&&... args)
      -> std::future<std::invoke_result_t<F, Args...>> {
    using R = std::invoke_result_t<F, Args...>;

    auto task_ptr = std::make_shared<std::packaged_task<R()>>(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...));

    enqueue_([task_ptr] { (*task_ptr)(); });
    return task_ptr->get_future();
  }

  template <typename F>
  void parallel_region(std::size_t n, F&& fn) {
    if (n == 0)
      return;

    // Avoid nested parallelism deadlocks/oversubscription.
    if (in_worker_) {
      fn(std::size_t{0});
      return;
    }

    const std::size_t tc = thread_count();
    if (tc == 0) {
      fn(std::size_t{0});
      return;
    }
    if (n > tc)
      n = tc;

    using Fn = std::decay_t<F>;
    struct Shared {
      std::atomic<std::size_t> remaining;
      std::mutex m;
      std::condition_variable cv;
      Fn fn;
      Shared(std::size_t n_, Fn&& f_) : remaining(n_), fn(std::move(f_)) {}
    };

    auto shared = std::make_shared<Shared>(n, Fn(std::forward<F>(fn)));

    for (std::size_t tid = 0; tid < n; ++tid) {
      enqueue_([shared, tid] {
        shared->fn(tid);
        if (shared->remaining.fetch_sub(1, std::memory_order_acq_rel) == 1) {
          std::lock_guard<std::mutex> lk(shared->m);
          shared->cv.notify_one();
        }
      });
    }

    std::unique_lock<std::mutex> lk(shared->m);
    shared->cv.wait(lk, [&] {
      return shared->remaining.load(std::memory_order_acquire) == 0;
    });
  }

 private:
  SimpleThreadPool() : done_(false) {
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0)
      hw = 2;
    const unsigned num_threads = hw;

    workers_.reserve(num_threads);
    for (unsigned i = 0; i < num_threads; ++i) {
      workers_.emplace_back([this] {
        // Per-worker AD tape (TLS) initialized once.
        static thread_local ChainableStack ad_tape;

        for (;;) {
          std::function<void()> task;
          {
            std::unique_lock<std::mutex> lock(mtx_);
            cv_.wait(lock, [&] { return done_ || !tasks_.empty(); });
            if (done_ && tasks_.empty())
              return;
            task = std::move(tasks_.front());
            tasks_.pop();
          }

          WorkerScope scope;  // sets in_worker_ for all tasks
          task();
        }
      });
    }
  }

  ~SimpleThreadPool() {
    {
      std::lock_guard<std::mutex> lock(mtx_);
      done_ = true;
    }
    cv_.notify_all();
    for (auto& th : workers_) {
      if (th.joinable())
        th.join();
    }
  }

  void enqueue_(std::function<void()> task) {
    {
      std::lock_guard<std::mutex> lock(mtx_);
      tasks_.emplace(std::move(task));
    }
    cv_.notify_one();
  }

  struct WorkerScope {
    WorkerScope() : prev_(in_worker_) { in_worker_ = true; }
    ~WorkerScope() { in_worker_ = prev_; }
    bool prev_;
  };

  static inline thread_local bool in_worker_ = false;

  std::vector<std::thread> workers_;
  std::queue<std::function<void()>> tasks_;
  std::mutex mtx_;
  std::condition_variable cv_;
  bool done_;
};

}  // namespace math
}  // namespace stan

#endif
