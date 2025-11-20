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
    
    //ThreadPool& pool = ThreadPool::instance();



    class WorkStealingThreadPool {
    public:
      // ----- Singleton access -----
      static WorkStealingThreadPool& instance() {
	static WorkStealingThreadPool inst;
	return inst;
      }

      WorkStealingThreadPool(const WorkStealingThreadPool&) = delete;
      WorkStealingThreadPool& operator=(const WorkStealingThreadPool&) = delete;

      // ----- Submit a task, get a future -----
      template <typename F, typename... Args>
      auto submit(F&& f, Args&&... args)
	-> std::future<std::invoke_result_t<F, Args...>> {
	using R = std::invoke_result_t<F, Args...>;

	auto task_ptr =
	  std::make_shared<std::packaged_task<R()>>(
						    std::bind(std::forward<F>(f), std::forward<Args>(args)...));

	std::function<void()> wrapper = [task_ptr]() { (*task_ptr)(); };

	std::size_t idx = worker_index;
	if (idx == kNoWorker) {
	  // Called from outside the pool: pick a worker round-robin.
	  idx = next_worker.fetch_add(1, std::memory_order_relaxed) %
            workers.size();
	}

	{
	  // Increase global task count *before* making it visible to workers
	  tasks_total.fetch_add(1, std::memory_order_relaxed);

	  std::lock_guard<std::mutex> lock(workers[idx]->mutex);
	  workers[idx]->tasks.emplace_back(std::move(wrapper));
	}

	// Wake one waiter (or do nothing if no one is waiting)
	cv.notify_one();

	return task_ptr->get_future();
      }

    private:
      struct Worker {
	std::deque<std::function<void()>> tasks;
	std::mutex mutex;
      };

      WorkStealingThreadPool()
	: done(false), next_worker(0), tasks_total(0) {
	std::size_t n = std::thread::hardware_concurrency();
	if (n == 0) {
	  n = 1;
	}

	workers.reserve(n);
	for (std::size_t i = 0; i < n; ++i) {
	  workers.emplace_back(std::make_unique<Worker>());
	}

	threads.reserve(n);
	for (std::size_t i = 0; i < n; ++i) {
	  threads.emplace_back([this, i] {
	    worker_index = i;
	    worker_loop(i);
	  });
	}
      }

      ~WorkStealingThreadPool() {
	done.store(true, std::memory_order_relaxed);
	cv.notify_all();

	for (auto& t : threads) {
	  if (t.joinable()) {
	    t.join();
	  }
	}
      }

      void worker_loop(std::size_t index) {
	while (true) {
	  std::function<void()> task;

	  if (pop_local_task(index, task) || steal_task(index, task)) {
	    // We successfully got a task, account for it and run it.
	    tasks_total.fetch_sub(1, std::memory_order_relaxed);
	    task();
	    continue;
	  }

	  // No work found right now. Go to sleep until work appears or shutdown.
	  std::unique_lock<std::mutex> lock(cv_mutex);
	  cv.wait(lock, [this] {
	    return done.load(std::memory_order_relaxed) ||
	      tasks_total.load(std::memory_order_relaxed) > 0;
	  });

	  if (done.load(std::memory_order_relaxed) &&
	      tasks_total.load(std::memory_order_relaxed) == 0) {
	    // Shutdown with no work remaining.
	    return;
	  }
	  // Otherwise: loop and try to grab work again.
	}
      }

      bool pop_local_task(std::size_t index, std::function<void()>& task) {
	Worker& w = *workers[index];
	std::lock_guard<std::mutex> lock(w.mutex);
	if (w.tasks.empty()) {
	  return false;
	}
	// Own worker: LIFO
	task = std::move(w.tasks.back());
	w.tasks.pop_back();
	return true;
      }

      bool steal_task(std::size_t thief_index, std::function<void()>& task) {
	const std::size_t n = workers.size();

	for (std::size_t i = 0; i < n; ++i) {
	  std::size_t victim_index = (thief_index + i + 1) % n;
	  Worker& w = *workers[victim_index];

	  std::lock_guard<std::mutex> lock(w.mutex);
	  if (!w.tasks.empty()) {
	    // Steal from front: FIFO
	    task = std::move(w.tasks.front());
	    w.tasks.pop_front();
	    return true;
	  }
	}
	return false;
      }

      std::vector<std::thread> threads;
      std::vector<std::unique_ptr<Worker>> workers;

      std::atomic<bool> done;
      std::atomic<std::size_t> next_worker;
      std::atomic<std::size_t> tasks_total;  // total tasks enqueued but not yet taken

      std::condition_variable cv;
      std::mutex cv_mutex;

      static constexpr std::size_t kNoWorker =
	std::numeric_limits<std::size_t>::max();

      static inline thread_local std::size_t worker_index = kNoWorker;
    };
  }  // namespace math
}  // namespace stan
#endif
