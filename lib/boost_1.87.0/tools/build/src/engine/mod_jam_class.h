/*
Copyright 2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_MOD_JAM_CLASS_H
#define B2_MOD_JAM_CLASS_H

#include "config.h"

#include "bind.h"
#include "lists.h"
#include "value.h"

#include <tuple>

/* tag::reference[]

[[b2.reference.modules.class]]
= `class` module.

end::reference[] */

namespace b2 { namespace jam { namespace klass {

/* tag::reference[]

== `b2::class::make`

====
[horizontal]
Jam:: `rule new ( class args * : * )`
{CPP}:: `std::string make(std::tuple<value_ref, list_ref> name_arg1, const lists
& rest, bind::context_ref_ context_ref)`
====

Instantiates a new instance of the given class and calls the `__init__` with
the given arguments. Returns the instance ID.

end::reference[] */
std::string make(std::tuple<value_ref, list_ref> name_arg1,
	const lists & rest,
	bind::context_ref_ context_ref);

/* tag::reference[]

== `b2::class::bases`

====
[horizontal]
Jam:: `rule bases ( class )`
{CPP}:: `list_ref bases(std::string class_name);`
====

Returns the base classes of the given class.

end::reference[] */
list_ref bases(std::string class_name);

/* tag::reference[]

== `b2::class::is_derived`

====
[horizontal]
Jam:: `rule is-derived ( class : bases + )`
{CPP}:: `bool is_derived(value_ref class_name, list_cref class_bases);`
====

Returns true when the given class is derived from the given bases.

end::reference[] */
bool is_derived(value_ref class_name, list_cref class_bases);

/* tag::reference[]

== `b2::class::is_instance`

====
[horizontal]
Jam:: `rule is-instance ( value )`
{CPP}:: `bool is_instance(std::string value);`
====

Returns true if the given value is an instance of any class.

end::reference[] */
bool is_instance(std::string value);

/* tag::reference[]

== `b2::class::is_a`

====
[horizontal]
Jam:: `rule is-a ( instance : type )`
{CPP}:: `bool is_a(std::string instance, value_ref type);`
====

Returns true if the given instance is of the given type.

end::reference[] */
bool is_a(std::string instance, value_ref type);

struct class_module : b2::bind::module_<class_module>
{
	const char * module_name = "class";
	static const char * init_code;

	template <class Binder>
	void def(Binder & binder)
	{
		binder.def(&make, "new", ("class" * _1 + "args" * _n) | "rest" * _r)
			.def(&bases, "bases", "class" * _1)
			.def(&is_derived, "is-derived", "class" * _1 | "bases" * _1n)
			.def(&is_instance, "is-instance", "value" * _1)
			.def(&is_a, "is-a", "instance" * _1 | "type" * _1);
		binder.eval(init_code);
		binder.loaded();
	}
};

}}} // namespace b2::jam::klass

#endif
