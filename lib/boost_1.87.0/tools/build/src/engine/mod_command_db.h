/*
Copyright 2024 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_COMMAND_DB_H
#define B2_COMMAND_DB_H

#include "config.h"

#include "bind.h"
#include "value.h"

namespace lyra {
class cli;
} // namespace lyra

namespace b2 {

namespace command_db {

void declare_args(lyra::cli &);

void set_output_dir(value_ref dirname);

} // namespace command_db

struct command_db_module : b2::bind::module_<command_db_module>
{
	const char * module_name = "command-db";

	template <class Binder>
	void def(Binder & binder)
	{
		binder.def(
			&command_db::set_output_dir, "set-output-dir", ("dirname" * _1));
		binder.loaded();
	}
};

} // namespace b2

#endif
