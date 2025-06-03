/*
Copyright 2022-2023 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "tasks.h"
#include "config.h"

#include "jam.h"
#include "mod_sysinfo.h"
#include "output.h"
#include "types.h"

#include <queue>
#include <vector>
#if B2_USE_STD_THREADS
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
#endif

namespace b2 { namespace task {

struct sync
{
	inline sync();
	inline void wait();
	inline bool signal();

#if B2_USE_STD_THREADS
	std::mutex arrived_mx;
	std::condition_variable arrived_cv;
	std::atomic_bool wait_arrived;
	std::atomic_bool signal_arrived;
#endif
};

inline sync::sync()
{
#if B2_USE_STD_THREADS
	wait_arrived = false;
	signal_arrived = false;
#endif
}

inline void sync::wait()
{
#if B2_USE_STD_THREADS
	// Indicate that we waiting.
	wait_arrived = true;
	// Wait for the signal that we can proceed.
	std::unique_lock<std::mutex> lock(arrived_mx);
	arrived_cv.wait(lock, [this]() { return signal_arrived.load(); });
#endif
}

inline bool sync::signal()
{
#if B2_USE_STD_THREADS
	// Wait for wait() to get called.
	if (!wait_arrived.load()) return false;
	// Tell the waiter that we arrived.
	std::unique_lock<std::mutex> lock(arrived_mx);
	signal_arrived = true;
	arrived_cv.notify_one();
#endif
	return true;
}

/*
A group of tasks that run in parallel within a limit of parallelism. The
parallelism limit is enforced by only dequeuing calls when possible.
*/
struct group::implementation
{
	executor & exec;
	unsigned parallelism = 0;
	unsigned running = 0;
	std::queue<std::function<void()>> pending;
	unsigned stat_max_running = 0;
	unsigned stat_total = 0;
	unsigned stat_max_pending = 0;
	mutex_t mx;
	sync finished;

	inline implementation(executor & e, unsigned p)
		: exec(e)
		, parallelism(p)
	{}

	// Add a task call to run.
	inline void call_queue(std::function<void()> f);

	// Remove and return a call ro run. If no call is available or possible
	// and empty function is returned.
	inline std::function<void()> call_dequeue();

	// Wait for all queued task calls to complete.
	inline void wait();
};

/*
A global executor of tasks from a collection of task groups. There is an overall
parallelism limit that follow the general Jam (-j) limit by default. Tasks from
groups are run within this limit regardless of what the group limits are. That
is, it will run the minimum of the execution limit and the group limits.
*/
struct executor::implementation
{
	std::vector<std::shared_ptr<group>> groups;
	unsigned call_count = 0;
	unsigned running_count = 0;
	mutex_t mx;

#if B2_USE_STD_THREADS
	std::vector<std::thread> runners;
	std::condition_variable call_cv;
#endif // B2_USE_STD_THREADS

	// Construction starts the task calling threads within the parallelism
	// limit.
	inline implementation(unsigned parallelism);

	// Add a group of tasks to run.
	inline void push_group(std::shared_ptr<group> g);

	// Signal that new tasks are available, and hence that the task threads
	// should run them.
	inline void call_signal();

	// Get a task call to execute from the collective set of group tasks.
	inline std::function<void()> call_get();

	// Indicate that an task call completed.
	inline void call_done();

	// Should we keep running (true), or complete ASAP (false).
	inline bool is_running();

	// Stop the threads in the execution pool regardless of pending tasks.
	inline void stop();

	// Waits for task calls to execute. This is teh body of each thread.
	void runner();
};

inline void group::implementation::call_queue(std::function<void()> f)
{
	// If we don't have parallel allotment. We opt to execute the call inline.
	if (parallelism == 0)
	{
		stat_total += 1;
		f();
		return;
	}
	// We wrap the bare task function to track how many are running at any time.
	// The running count is used in the wait and to enforce the group
	// parallelism limit.
	auto the_call = [this, f]() {
		{
			scope_lock_t lock(mx);
			running += 1;
			stat_max_running = std::max(running, stat_max_running);
			stat_total += 1;
		}
		f();
		bool signal_finished = false;
		{
			scope_lock_t lock(mx);
			running -= 1;
			signal_finished = pending.empty() && running == 0;
		}
		if (signal_finished) finished.signal();
	};
	{
		scope_lock_t lock(mx);
		pending.push(std::move(the_call));
		stat_max_pending = std::max(unsigned(pending.size()), stat_max_pending);
	}
	// Signal the executor that there are tasks to run.
	exec.i->call_signal();
}

inline std::function<void()> group::implementation::call_dequeue()
{
	bool signal_finished = false;
	std::function<void()> result;
	{
		scope_lock_t lock(mx);
		signal_finished = pending.empty() && running == 0;
		// We only return tasks when we have them, and when we have enough
		// parallelism.
		if (!pending.empty() && running < parallelism)
		{
			result = std::move(pending.front());
			pending.pop();
		}
	}
	if (signal_finished) finished.signal();
	return result;
}

inline void group::implementation::wait()
{
	// Without parallelism there's nothing to wait. As everything was executed
	// inline already.
	if (parallelism == 0) return;
	// Signal the tasks that we are waiting for completion.
	exec.i->call_signal();
	// Wait for completion of the group tasks.
	finished.wait();
}

inline executor::implementation::implementation(unsigned parallelism)
{
	// No need to launch anything if we aren't parallel.
	if (parallelism == 0) return;
	// Launch the threads to cover the expected parallelism.
	scope_lock_t lock(mx);
#if B2_USE_STD_THREADS
	runners.reserve(parallelism);
	for (; parallelism > 0; --parallelism)
	{
		try
		{
			runners.emplace_back([this]() { runner(); });
			running_count += 1;
		}
		catch (const std::system_error & e)
		{
			err_printf("Task startup error: %s", e.what());
		}
	}
#endif
}

inline void executor::implementation::push_group(std::shared_ptr<group> g)
{
	scope_lock_t lock(mx);
	groups.push_back(g);
}

inline void executor::implementation::call_signal()
{
	scope_lock_t lock(mx);
#if B2_USE_STD_THREADS
	call_cv.notify_one();
#endif
}

inline std::function<void()> executor::implementation::call_get()
{
	std::function<void()> result;
	scope_lock_t lock(mx);
#if B2_USE_STD_THREADS
	// We only dequeue task calls when we have a thread to run them.
	if (call_count < runners.size())
	{
		for (auto & group : groups)
		{
			result = group->i->call_dequeue();
			if (result)
			{
				call_count += 1;
				return result;
			}
		}
	}
	// We don't have tasks to run, wait for some to become available.
	call_cv.wait(lock);
#endif
	return result;
}

inline void executor::implementation::call_done()
{
	scope_lock_t lock(mx);
	call_count -= 1;
}

inline bool executor::implementation::is_running()
{
	scope_lock_t lock(mx);
	return running_count > 0;
}

inline void executor::implementation::stop()
{
#if B2_USE_STD_THREADS
	// Stop all the runner threads (i.e. signal and wait).
	std::vector<std::thread> to_join;
	{
		scope_lock_t lock(mx);
		running_count = 0;
		to_join.swap(runners);
	}
	for (auto & t : to_join)
	{
		{
			scope_lock_t lock(mx);
			call_cv.notify_all();
		}
		t.join();
	}
#endif
}

void executor::implementation::runner()
{
	while (is_running())
	{
		auto f = call_get();
		if (f)
		{
			try
			{
				f();
			}
			catch (const std::exception & e)
			{
				err_printf("Task runner function error: %s", e.what());
			}
			call_done();
		}
	}
}

namespace {
unsigned get_parallelism(int parallelism)
{
#if B2_USE_STD_THREADS
	return parallelism >= 0
		? parallelism
		: std::min(unsigned(globs.jobs), system_info().cpu_thread_count());
#else
	return 0;
#endif
}
} // namespace

group::group(executor & exec, int parallelism)
	: i(std::make_shared<implementation>(exec, parallelism))
{}

group::~group()
{
	// out_printf("b2::task::group.. MAX = %u, TOTAL = %u, MAX PENDING = %u\n",
	// 	i->stat_max_running, i->stat_total, i->stat_max_pending);
}

void group::queue(std::function<void()> && f) { i->call_queue(f); }

void group::wait() { i->wait(); }

executor::executor(int parallelism)
	: i(std::make_shared<implementation>(get_parallelism(parallelism)))
{}

executor::~executor() { i->stop(); }

std::shared_ptr<group> executor::make(int parallelism)
{
	auto result = std::make_shared<group>(*this,
		std::min(unsigned(i->runners.size()), get_parallelism(parallelism)));
	i->push_group(result);
	return result;
}

executor & executor::get()
{
	static executor e;
	return e;
}

}} // namespace b2::task
