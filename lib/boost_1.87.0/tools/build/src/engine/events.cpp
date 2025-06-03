/*
Copyright 2024 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "events.h"

#include <algorithm>
#include <memory>
#include <set>
#include <vector>

namespace b2 {

struct event_base
{
	event_tag tag;
	int32_t priority;
	uint64_t id;

	inline bool operator<(const event_base & o) const
	{
		if (tag < o.tag) return true;
		if ((0 - priority) < (0 - o.priority)) return true;
		return id < o.id;
	}
};

template <typename F>
struct event : public event_base
{
	F call;
};

struct events
{
	struct ecmp
	{
		inline bool operator()(const event_base * a, const event_base * b) const
		{
			return (*a) < (*b);
		}
	};
	std::set<const event_base *, ecmp> sorted_items;
	std::vector<std::unique_ptr<event_base>> items;

	static events & get()
	{
		static events e;
		return e;
	}

	template <typename F>
	uint64_t add(event_tag tag, F && call, int32_t priority)
	{
		uint64_t id = items.empty() ? 1 : items.back()->id + 1;
		std::unique_ptr<event<F>> e(new event<F>);
		e->tag = tag;
		e->priority = priority;
		e->id = id;
		e->call = call;
		sorted_items.insert(e.get());
		items.push_back(std::move(e));
		return id;
	}

	void remove(uint64_t e)
	{
		auto i = std::lower_bound(items.begin(), items.end(), e,
			[](const std::unique_ptr<event_base> & a, uint64_t id) -> bool {
				return a->id < id;
			});
		if (i != items.end() && (*i)->id == e)
		{
			items.erase(i);
		}
	}

	template <typename... Args>
	void trigger(event_tag tag, Args... args)
	{
		using E = event<std::function<void(Args...)>>;
		static event_base x = { tag, std::numeric_limits<int32_t>::max(), 0 };
		auto i = sorted_items.lower_bound(&x);
		static event_base y = { tag, std::numeric_limits<int32_t>::min(), 0 };
		auto j = sorted_items.lower_bound(&y);
		if (j != sorted_items.end()) ++j;
		for (; i != j; ++i)
		{
			if ((*i)->tag != tag) break;
			static_cast<const E *>(*i)->call(args...);
		}
	}
};

void remove_event_callback(uint64_t e) { events::get().remove(e); }

template <>
uint64_t add_event_callback(
	event_tag tag, std::function<void(TARGET *)> && call, int32_t priority)
{
	return events::get().add(tag, std::move(call), priority);
}

void trigger_event_pre_exec_cmd(TARGET * t)
{
	events::get().trigger<TARGET *>(event_tag::pre_exec_cmd, t);
}

template <>
uint64_t add_event_callback(
	event_tag tag, std::function<void(int)> && call, int32_t priority)
{
	return events::get().add(tag, std::move(call), priority);
}

void trigger_event_exit_main(int status)
{
	events::get().trigger<int>(event_tag::exit_main, status);
}

} // namespace b2
