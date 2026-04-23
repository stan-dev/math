/*
Copyright 2019-2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_MOD_STRING_H
#define B2_MOD_STRING_H

#include "config.h"

#include "bind.h"
#include "lists.h"
#include "optval.h"
#include "value.h"

#include <string>
#include <vector>

/* tag::reference[]

[[b2.reference.modules.string]]
= `string` module.

end::reference[] */

namespace b2 {
/* tag::reference[]

== `b2::string_whitespace`

====
[horizontal]
Jam:: `rule whitespace ( )`
{CPP}:: `b2::value_ref string_whitespace();`
====

Returns the canonical set of whitespace characters, as a single string.

end::reference[] */
value_ref string_whitespace();

/* tag::reference[]

== `b2::string_chars`

====
[horizontal]
Jam:: `rule chars ( string )`
{CPP}:: `b2::list_ref string_chars(b2::value_ref s);`
====

Splits the given `string` into a list of strings composed of each character of
the `string` in sequence.

end::reference[] */
list_ref string_chars(value_ref s);

/* tag::reference[]

== `b2::string_abbreviate`

====
[horizontal]
Jam:: `rule abbreviate ( string )`
{CPP}:: `b2::value_ref string_abbreviate(b2::value_ref s);`
====

Apply a set of standard transformations to string to produce an abbreviation
no more than 5 characters long.

end::reference[] */
value_ref string_abbreviate(value_ref s);

/* tag::reference[]

== `b2::string_join`

====
[horizontal]
Jam:: `rule join ( strings * : separator ? )`
{CPP}:: `b2::value_ref string_join(b2::list_cref strings, b2::value_ref
separator);`
====

Concatenates the given `strings`, inserting the given `separator` between each
string.

end::reference[] */
value_ref string_join(list_cref strings, value_ref separator);

/* tag::reference[]

== `b2::string_words`

====
[horizontal]
Jam:: `rule words ( string : whitespace * )`
{CPP}:: `b2::list_ref string_words(std::string s, b2::list_cref whitespace);`
====

Split a string into whitespace separated words.

end::reference[] */
list_ref string_words(std::string s, list_cref whitespace);

/* tag::reference[]

== `b2::string_is_whitespace`

====
[horizontal]
Jam:: `rule is-whitespace ( string ? )`
{CPP}:: `bool string_is_whitespace(b2::value_ref s);`
====

Check that the given `string` is composed entirely of whitespace.

end::reference[] */
bool string_is_whitespace(value_ref s);

struct string_module : b2::bind::module_<string_module>
{
	const char * module_name = "string";
	static const char * init_code;

	template <class Binder>
	void def(Binder & binder)
	{
		binder.def(&b2::string_whitespace, "whitespace")
			.def(&b2::string_chars, "chars", "string" * _1)
			.def(&b2::string_abbreviate, "abbreviate", "string" * _1)
			.def(&b2::string_join, "join", "strings" * _n | "separator" * _01)
			.def(&b2::string_words, "words", "string" * _1 | "whitespace" * _n)
			.def(&b2::string_is_whitespace, "is-whitespace", "string" * _01);
		binder.eval(init_code);
		binder.loaded();
	}
};
} // namespace b2

#endif
