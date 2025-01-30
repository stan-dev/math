/*
Copyright 2019-2023 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_MOD_REGEX_H
#define B2_MOD_REGEX_H

#include "config.h"

#include "bind.h"
#include "lists.h"
#include "optval.h"
#include "types.h"
#include "value.h"

#include <tuple>
#include <vector>

/* tag::reference[]

[[b2.reference.modules.regex]]
= `regex` module.

Contains rules for string processing using regular expressions.

* `"x*"` matches the pattern `"x"` zero or more times.
* `"x+"` matches `"x"` one or more times.
* `"x?"` matches `"x"` zero or one time.
* `"[abcd]"` matches any of the characters, `"a"`, `"b"`, `"c"`, and
`"d"`. A character range such as `"[a-z]"` matches any character between
`"a"` and `"z"`. `"[^abc]"` matches any character which is not `"a"`,
`"b"`, or `"c"`.
* `"x|y"` matches either pattern `"x"` or pattern `"y"`
* `(x)` matches `"x"` and captures it.
* `"^"` matches the beginning of the string.
* `"$"` matches the end of the string.
* `"<"` matches the beginning of a word.
* `">"` matches the end of a word.

See also: link:#jam.language.rules.builtins.utility.\_match__[MATCH]

end::reference[] */

namespace b2 {

/* tag::reference[]

== `b2::regex_split`

====
[horizontal]
Jam:: `rule split ( string separator )`
{CPP}:: `b2::list_ref rb2::egex_split(const std::tuple<b2::value_ref,
b2::value_ref> & string_separator);`
====

Returns a list of the following substrings:

. from beginning till the first occurrence of 'separator' or till the end,
. between each occurrence of 'separator' and the next occurrence,
. from the last occurrence of 'separator' till the end.

If no separator is present, the result will contain only one element.

end::reference[] */
list_ref regex_split(const std::tuple<value_ref, value_ref> & string_separator);

/* tag::reference[]

== `b2::regex_split_each`

====
[horizontal]
Jam:: `rule split-list ( list * : separator )`
{CPP}:: `b2::list_ref b2::regex_split_each(b2::list_cref to_split, b2::value_ref
separator);`
====

Returns the concatenated results of Applying regex.split to every element of
the list using the separator pattern.

end::reference[] */
list_ref regex_split_each(list_cref to_split, value_ref separator);

/* tag::reference[]

== `b2::regex_match`

====
[horizontal]
Jam:: `rule match ( pattern : string : indices * )`
{CPP}:: `b2::list_ref regex_match(b2::value_ref pattern, b2::value_ref string,
const std::vector<int_t> & indices);`
====

Match `string` against `pattern`, and return the elements indicated by
`indices`.

end::reference[] */
list_ref regex_match(
	value_ref pattern, value_ref string, const std::vector<int_t> & indices);

/* tag::reference[]

== `b2::regex_transform`

====
[horizontal]
Jam:: `rule transform ( list * : pattern : indices * )`
{CPP}:: `b2::list_ref regex_transform(b2::list_cref list, b2::value_ref pattern,
const std::vector<int_t> & indices);`
====

Matches all elements of `list` against the `pattern` and returns a list of
elements indicated by indices of all successful matches. If `indices` is
omitted returns a list of first parenthesized groups of all successful
matches.

end::reference[] */
list_ref regex_transform(
	list_cref list, value_ref pattern, const std::vector<int_t> & indices);

/* tag::reference[]

== `b2::regex_escape`

====
[horizontal]
Jam:: `rule escape ( string : symbols : escape-symbol )`
{CPP}:: `b2::value_ref regex_escape(b2::value_ref string,b2:: value_ref symbols,
b2::value_ref escape_symbol);`
====

Escapes all of the characters in `symbols` using the escape symbol
`escape-symbol` for the given string, and returns the escaped string.

end::reference[] */
value_ref regex_escape(
	value_ref string, value_ref symbols, value_ref escape_symbol);

/* tag::reference[]

== `b2::regex_replace`

====
[horizontal]
Jam:: `rule replace ( string match replacement )`
{CPP}:: `b2::value_ref regex_replace(const std::tuple<b2::value_ref,
b2::value_ref, b2::value_ref> & string_match_replacement);`
====

Replaces occurrences of a match string in a given string and returns the new
string. The match string can be a regex expression.

* `string` -- The string to modify.
* `match` -- The characters to replace.
* `replacement` -- The string to replace with.

end::reference[] */
value_ref regex_replace(const std::tuple<value_ref, value_ref, value_ref> &
		string_match_replacement);

/* tag::reference[]

== `b2::regex_replace_each`

====
[horizontal]
Jam:: `rule replace-list ( list * : match : replacement )`
{CPP}:: `b2::list_ref regex_replace_each(b2::list_cref list, b2::value_ref
match, b2::value_ref replacement);`
====

Replaces occurrences of a match string in a given list of strings and returns
a list of new strings. The match string can be a regex expression.

* `list` -- The list of strings to modify.
* `match` -- The search expression.
* `replacement` -- The string to replace with.

end::reference[] */
list_ref regex_replace_each(
	list_cref list, value_ref match, value_ref replacement);

/* tag::reference[]

== `b2::regex_grep`

====
[horizontal]
Jam:: `rule grep ( directories + : files + : patterns + : result_expressions *
: options * )`
{CPP}:: `b2::list_ref regex_grep(b2::list_cref directories, b2::list_cref
files, b2::list_cref patterns, list_cref result_expressions,
list_cref options);`
====

Match any of the `patterns` against the globbed `files` in `directories`, and
return a list of files and indicated `result_expressions` (file1, re1, re..,
...). The `result_expressions` are indices from `0` to `10`. Where `0` is the
full match.

end::reference[] */
list_ref regex_grep(list_cref directories,
	list_cref files,
	list_cref patterns,
	list_cref result_expressions,
	list_cref options);

struct regex_module : b2::bind::module_<regex_module>
{
	const char * module_name = "regex";
	static const char * init_code;

	template <typename Binder>
	void def(Binder & binder)
	{
		binder.def(&regex_split, "split", "string" * _1 + "separator" * _1)
			.def(
				&regex_split_each, "split-list", "list" * _n | "separator" * _1)
			.def(&regex_match, "match",
				"pattern" * _1 | "string" * _1 | "indices" * _n)
			.def(&regex_transform, "transform",
				"list" * _n | "pattern" * _1 | "indices" * _n)
			.def(&regex_escape, "escape",
				"string" * _1 | "symbols" * _1 | "escape-symbol" * _1)
			.def(&regex_replace, "replace",
				"string" * _1 + "match" * _1 + "replacement" * _1)
			.def(&regex_replace_each, "replace-list",
				"list" * _n | "match" * _1 | "replacement" * _1)
			.def(&regex_grep, "grep",
				"directories" * _1n | "files" * _1n | "patterns" * _1n
					| "result_expressions" * _n | "options" * _n);
		binder.eval(init_code);
		binder.loaded();
	}
};
} // namespace b2

#endif
