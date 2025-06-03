/*
Copyright 2022-2023 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_TASKS_H
#define B2_TASKS_H

#include "config.h"

#include "mod_sysinfo.h"

#include <functional>
#include <memory>

/*
= Tasks

Utility classes for parallel invocation of "tasks" (i.e. functions).
*/

namespace b2 { namespace task {

class executor;

/*
== `b2::task::group`

A task group is a collection of tasks that will execute within a parallelism
limit of their own and can be waited on collectively to complete.
*/
class group
{
	public:
	group(executor & exec, int parallelism = -1);
	~group();

	// Add a task function to execute async. The functions are executed in no
	// particular order.
	void queue(std::function<void()> && f);

	// Wait for all the task functions in the group to complete.
	void wait();

	private:
	struct implementation;
	std::shared_ptr<implementation> i;

	friend class executor;
};

/*
== `b2::task::executor`

Global task execution queue that has a global parallelism limit. By default the
parallelism limit matches the `b2::system_info::cpu_thread_count` value.
*/
class executor
{
	public:
	executor(int parallelism = -1);
	~executor();

	// The global executor instance.
	static executor & get();

	// Create a task group.
	std::shared_ptr<group> make(int parallelism = -1);

	private:
	struct implementation;
	std::shared_ptr<implementation> i;

	friend class group;
};

}} // namespace b2::task

#endif
