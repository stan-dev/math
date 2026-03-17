/*
Copyright 2019-2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_MOD_SET_H
#define B2_MOD_SET_H

#include "config.h"

#include "bind.h"
#include "lists.h"
#include "optval.h"
#include "value.h"

#include <string>
#include <unordered_set>
#include <vector>

/* tag::reference[]

= `set` module.

Classes and functions to manipulate sets of unique values.

== `b2::set`

Set class contains unique values.

=== `b2::set::add`

====
[horizontal]
Jam:: `rule add ( elements * )`
{CPP}:: `void b2::set::add(b2::list_cref elements);`, `void b2::set::add(const
b2::set & elements);`
====

Add the `elements` to the set.

=== `b2::set::contains`

====
[horizontal]
Jam:: `rule contains ( element )`
{CPP}:: `bool b2::set::contains(b2::value_ref element) const;`
====

Does the set contain the given `element`.

=== `b2::set::to_list`

====
[horizontal]
Jam:: `rule list ( )`
{CPP}:: `b2::list_ref b2::set::to_list() const;`
====

Return a list with all the elements of the set.

=== `b2::set::difference`

====
[horizontal]
Jam:: `rule difference ( set1 * : set2 * )`
{CPP}:: `static b2::list_ref b2::set::difference(b2::list_cref set1,
b2::list_cref set2);`
====

Returns the elements of `set1` that are not in `set2`.

=== `b2::set::intersection`

====
[horizontal]
Jam:: `rule intersection ( set1 * : set2 * )`
{CPP}:: `static b2::list_ref b2::set::intersection(b2::list_cref set1,
b2::list_cref set2);`
====

Removes all the items appearing in both `set1` & `set2`.

=== `b2::set::equal`

====
[horizontal]
Jam:: `rule equal ( set1 * : set2 * )`
{CPP}:: `static bool b2::set::equal(b2::list_cref set1, b2::list_cref set2);`
====

Returns whether `set1` & `set2` contain the same elements. Note that this
ignores any element ordering differences as well as any element
duplication.

end::reference[] */

namespace b2 {

class set : public object
{
	public:
	void add(list_cref values);
	void add(const set & value);
	bool contains(value_ref value) const;
	list_ref to_list() const;
	static list_ref difference(list_cref a, list_cref b);
	static list_ref intersection(list_cref a, list_cref b);
	static bool equal(list_cref a, list_cref b);

	private:
	using value_set = std::unordered_set<value_ref,
		value_ref::hash_function,
		value_ref::equal_function>;

	value_set elements;
};

struct set_module : b2::bind::module_<set_module>
{
	const char * module_name = "set";
	static const char * init_code;

	template <class Binder>
	void def(Binder & binder)
	{
		binder.def(&set::difference, "difference", "set1" * _n | "set2" * _n)
			.def(&set::intersection, "intersection", "set1" * _n | "set2" * _n)
			.def(&set::equal, "equal", "set1" * _n | "set2" * _n);
		binder.def_class("set", type_<set>())
			.def(init_<>())
			.def(static_cast<void (set::*)(b2::list_cref)>(&set::add), "add",
				"elemets" * _n)
			.def(&set::contains, "contains", "element" * _1)
			.def(&set::to_list, "list");
		binder.eval(init_code);
		binder.loaded();
	}
};

} // namespace b2

#endif
