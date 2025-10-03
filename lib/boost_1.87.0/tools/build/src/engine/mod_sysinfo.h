/*
Copyright 2019-2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_MOD_SYSINFO_H
#define B2_MOD_SYSINFO_H

#include "config.h"

#include "bind.h"
#include "value.h"

#include <string>

namespace b2 {

/*
Provides information about the system, hardware and software, we are
running in.
*/
class system_info : public object
{
	public:
	system_info();

	/*
	Returns the number of physical CPU cores if available. Otherwise
	returns 1.

	Currently implemented for: OS_MACOSX.
	*/
	unsigned int cpu_core_count();

	/*
	Returns the number of logical CPU threads is available. Otherwise
	returns `spu_core_count()`.

	Currently implemented for: OS_MACOSX.
	*/
	unsigned int cpu_thread_count();

	private:
	unsigned int cpu_core_count_ = 0;
	unsigned int cpu_thread_count_ = 0;
};

struct stacktrace
{
	static std::string to_string();
};

struct sysinfo_module : b2::bind::module_<sysinfo_module>
{
	const char * module_name = "sysinfo";

	template <class Binder>
	void def(Binder & binder)
	{
		binder.def_class("system_info", type_<b2::system_info>())
			.def(init_<>())
			.def(&b2::system_info::cpu_core_count, "cpu_core_count")
			.def(&b2::system_info::cpu_thread_count, "cpu_thread_count");
		binder.loaded();
	}
};

} // namespace b2

#endif
