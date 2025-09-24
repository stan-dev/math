/*
Copyright 2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_MOD_PATH_H
#define B2_MOD_PATH_H

#include "config.h"

#include "bind.h"
#include "lists.h"
#include "value.h"

/* tag::reference[]

[[b2.reference.modules.path]]
= `path` module.

end::reference[] */

namespace b2 { namespace paths {

/* tag::reference[]

== `b2::paths::exists`

====
[horizontal]
Jam:: `rule exists ( file )`
{CPP}:: `bool exists(value_ref location);`
====

Returns true is the specified file exists.

end::reference[] */
bool exists(value_ref location);

/* tag::reference[]

== `b2::paths::normalize_all`

====
[horizontal]
Jam:: `rule normalize-all ( paths * )`
{CPP}:: `list_ref normalize_all(list_cref paths);`
====

Transform each path in the list, with all backslashes converted to forward
slashes and all detectable redundancy removed.

end::reference[] */
list_ref normalize_all(list_cref paths);

struct paths_module : b2::bind::module_<paths_module>
{
	const char * module_name = "path";

	template <class Binder>
	void def(Binder & binder)
	{
		binder.def(&exists, "exists", "location" * _1)
			.def(&normalize_all, "normalize-all", "paths" * _n);
	}
};

}} // namespace b2::paths

#endif
