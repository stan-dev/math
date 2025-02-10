/*
Copyright 2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_MOD_BUILTIN_H
#define B2_MOD_BUILTIN_H

#include "config.h"

#include "bind.h"
#include "value.h"

namespace b2 { namespace jam { namespace builtin {

/* tag::reference[]

== `b2::require_b2`

====
[horizontal]
Jam:: `rule require-b2 ( minimum : maximum ? )`
{CPP}:: `void require_b2(value_ref minimum, value_ref maximum,
bind::context_ref_ context_ref);`
====

Checks that the b2 engine version is at least `minimum` and strictly less than
`maximum`. If no `maximum` is given the version is matched to be at least
the `minimum`. If the check fails b2 will exit with an error message and
failure.

end::reference[] */
void require_b2(
	value_ref minimum, value_ref maximum, bind::context_ref_ context_ref);

struct root_module : b2::bind::module_<root_module>
{
	const char * module_name = nullptr;
	static const char * init_code;

	template <class Binder>
	void def(Binder & binder)
	{
		binder.def(&require_b2, "require-b2", "minimum" * _1 | "maximum" * _01);
		binder.eval(init_code);
	}
};

}}} // namespace b2::jam::builtin

#endif
