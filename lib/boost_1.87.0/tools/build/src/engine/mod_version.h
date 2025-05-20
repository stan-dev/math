/*
Copyright 2019-2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_MOD_VERSION_H
#define B2_MOD_VERSION_H

#include "bind.h"
#include "config.h"
#include <vector>

namespace b2 {
/* tag::reference[]

[[b2.reference.modules.version]]
= `version` module.

end::reference[] */

/* tag::reference[]

== `b2::version_less`

====
[horizontal]
Jam:: `rule version-less ( lhs + : rhs + )`
{CPP}:: `bool version_less(const std::vector<int> & lhs, const std::vector<int>
& rhs);`
====

Returns `true` if the first version, `lhs`, is semantically less than the
second version, `rhs`.

end::reference[] */
bool version_less(const std::vector<int> & lhs, const std::vector<int> & rhs);

struct version_module : b2::bind::module_<version_module>
{
	const char * module_name = "version";

	template <class Binder>
	void def(Binder & binder)
	{
		binder.def(
			&b2::version_less, "version-less", "lhs" * _1n | "rhs" * _1n);
	}
};
} // namespace b2

#endif
