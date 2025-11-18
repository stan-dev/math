#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>

#include <vector>
#include <thread>
#include <future>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>



namespace stan {
  namespace math {
    class chainable_alloc;
    class vari_base;
    using ChainableStack = AutodiffStackSingleton<vari_base, chainable_alloc>;

    class ThreadPool {
    public:
      // ----- Access to the global thread pool -----
      static ThreadPool& instance() {
        static ThreadPool inst;   // created on first call, thread-safe
        return inst;
      }

      // ----- Submit a task and get a future -----
      template <class F, class... Args>
      auto submit(F&& f, Args&&... args)
        -> std::future<std::invoke_result_t<F, Args...>>
      {
        using R = std::invoke_result_t<F, Args...>;

        auto task = std::make_shared<std::packaged_task<R()>>(
							      std::bind(std::forward<F>(f), std::forward<Args>(args)...));

        {
	  std::unique_lock<std::mutex> lock(mtx);
	  tasks.emplace([task]{ (*task)(); });
        }

        cv.notify_one();
        return task->get_future();
      }

      // ----- Prevent copying -----
      ThreadPool(const ThreadPool&) = delete;
      ThreadPool& operator=(const ThreadPool&) = delete;

    private:
      ThreadPool() : stop(false) {
        const size_t n = std::thread::hardware_concurrency();
        for (size_t i = 0; i < n; ++i) {
	  workers.emplace_back([this] {
	    for (;;) {
	      std::function<void()> task;

	      {
		std::unique_lock<std::mutex> lock(mtx);
		cv.wait(lock, [&]{
		  return stop || !tasks.empty();
		});
		if (stop && tasks.empty())
		  return;
		task = std::move(tasks.front());
		tasks.pop();
	      }

	      task();
	    }
	  });
        }
      }

      ~ThreadPool() {
        {
	  std::unique_lock<std::mutex> lock(mtx);
	  stop = true;
        }
        cv.notify_all();
        for (auto& t : workers)
	  t.join();
      }

      std::vector<std::thread> workers;
      std::queue<std::function<void()>> tasks;

      std::mutex mtx;
      std::condition_variable cv;
      bool stop;
    };
    
    ThreadPool& pool = ThreadPool::instance();
  
  }  // namespace math
}  // namespace stan
#endif
