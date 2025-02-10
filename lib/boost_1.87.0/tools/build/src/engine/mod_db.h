/*
Copyright 2024 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_MOD_DB_H
#define B2_MOD_DB_H

#include "config.h"

#include "bind.h"
#include "lists.h"
#include "value.h"

#include <map>

/* tag::reference[]

= `db` module.

Classes and functions to manage structured data.

end::reference[] */

namespace b2 {

/* tag::reference[]

== `b2::property_db` (`property-db`)

Container for values structured as a tree with keys that are the tree path to
the value. Arrays and objects (named fields) are supported.

end::reference[] */
struct property_db : public object
{
	void emplace(list_cref k, value_ref v);
	bool write_file(value_ref filename, value_ref format);
	std::string dump(value_ref format);

	private:
	value_ref array_tag { "[]" };
	struct key_less
	{
		bool operator()(const list_ref & a, const list_ref & b) const
		{
			for (list_ref::size_type i = 0; i < a.size() && i < b.size(); ++i)
			{
				value_ref a_v = a[i];
				value_ref b_v = b[i];
				int c = a_v->compare(*b_v);
				if (c < 0) return true;
				if (c > 0) return false;
			}
			return a.size() < b.size();
		}
	};
	using db_type = std::map<list_ref, value_ref, key_less>;
	db_type db;
	bool write_file_json(value_ref filename);
	std::string dump_json();
};

/* tag::reference[]

=== `b2::property_db::emplace`

====
[horizontal]
Jam:: `rule emplace ( key + : value )`
{CPP}:: `void emplace(list_cref k, value_ref v);`
====

Set, or add, element at path `key` with `value`. The path can contain two kinds
of position items: an array index or an object member. An array index is
`"[]" _n_`. Where the `"[]"` indicates that a zero based index in the array,
`_n_` follows. Anything else is treated as a member field name.

end::reference[] */

/* tag::reference[]

=== `b2::property_db::write_file`

====
[horizontal]
Jam:: `rule write-file ( filename : format ? )`
{CPP}:: `void write_file(value_ref filename, value_ref format);`
====

Writes out a representation of the data to the given `filename` file formatted
as `format`. Supported formats are: `JSON`.

end::reference[] */

/* tag::reference[]

=== `b2::property_db::dump`

====
[horizontal]
Jam:: `rule dump ( format ? )`
{CPP}:: `std::string dump(value_ref format);`
====

Writes out a representation of the data as a string formatted as `format`.
Supported formats are: `JSON`.

end::reference[] */

struct db_module : b2::bind::module_<db_module>
{
	const char * module_name = "db";
	static const char * init_code;

	template <class Binder>
	void def(Binder & binder)
	{
		binder.def_class("property-db", type_<property_db>())
			.def(init_<>())
			.def(&property_db::emplace, "emplace",
				("key" * _1n) | ("value" * _1))
			.def(&property_db::write_file, "write-file",
				("filename" * _1) | ("format" * _01))
			.def(&property_db::dump, "dump", ("format" * _01));
		binder.eval(init_code);
		binder.loaded();
	}
};

} // namespace b2

#endif
