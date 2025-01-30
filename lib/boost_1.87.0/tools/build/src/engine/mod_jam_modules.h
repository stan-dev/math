/*
Copyright 2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_MOD_MODULES_H
#define B2_MOD_MODULES_H

#include "config.h"

#include "bind.h"
#include "lists.h"
#include "value.h"

#include "mod_path.h"

#include <string>

/* tag::reference[]

[[b2.reference.modules.modules]]
= `modules` module.

The `modules` module defines basic functionality for handling modules.

A module defines a number of rules that can be used in other modules.
Modules can contain code at the top level to initialize the module. This
code is executed the first time the module is loaded.

NOTE: A Jamfile is a special kind of module which is managed by the build
system. Although they cannot be loaded directly by users, the other
features of modules are still useful for Jamfiles.

Each module has its own namespaces for variables and rules. If two
modules A and B both use a variable named X, each one gets its own copy
of X. They won't interfere with each other in any way. Similarly,
importing rules into one module has no effect on any other module.

Every module has two special variables. `$(__file__)` contains the name
of the file that the module was loaded from and `$(__name__)` contains
the name of the module.

NOTE: `$(__file__)` does not contain the full path to the file. If you need
this, use `modules.binding`.

end::reference[] */

namespace b2 { namespace jam { namespace modules {

/* tag::reference[]

== `b2::modules::binding`

====
[horizontal]
Jam:: `rule binding ( module_name )`
{CPP}:: `value_ref binding(std::string module_name);`
====

Returns the filesystem binding of the given module.

For example, a module can get its own location with:

[source,jam]
----
me = [ modules.binding $(__name__) ] ;
----

end::reference[] */
list_ref binding(std::string module_name);

/* tag::reference[]

== `b2::modules::record_binding`

====
[horizontal]
Jam:: `rule record-binding ( module_name : binding )`
{CPP}:: `void record_binding(std::string module_name, value_ref value);`
====

This helper is used by load to record the binding (path) of each loaded module.

WARNING: The `module_name` is ignored. Instead the internal tracking of the
currently loading module is used to record the binding.

end::reference[] */
void record_binding(std::string module_name, value_ref value);

/* tag::reference[]

== `b2::modules::poke`

====
[horizontal]
Jam:: `rule poke ( module_name ? : variables + : value * )`
{CPP}:: `void poke(std::string module_name, list_cref variables, list_cref
value);`
====

Sets the module-local value of a variable. This is the most reliable way to
set a module-local variable in a different module; it eliminates issues of
name shadowing due to dynamic scoping.

For example, to set a variable in the global module:

[source,jam]
----
modules.poke : ZLIB_INCLUDE : /usr/local/include ;
----

end::reference[] */
void poke(std::string module_name, list_cref variables, list_cref value);

/* tag::reference[]

== `b2::modules::peek`

====
[horizontal]
Jam:: `rule peek ( module_name ? : variables + )`
{CPP}:: `list_ref peek(std::string module_name, list_cref variables);`
====

Returns the module-local value of a variable. This is the most reliable way to
examine a module-local variable in a different module; it eliminates issues of
name shadowing due to dynamic scoping.

For example, to read a variable from the global module:

[source,jam]
----
local ZLIB_INCLUDE = [ modules.peek : ZLIB_INCLUDE ] ;
----

end::reference[] */
list_ref peek(std::string module_name, list_cref variables);

/* tag::reference[]

== `b2::modules::clone_rules`

====
[horizontal]
Jam:: `rule clone-rules ( source_module target_module )`
{CPP}:: `void clone_rules(std::tuple<std::string, std::string>
source_target_modules);`
====

Define exported copies in `target_module` of all rules exported from
`source_module`. Also make them available in the global module with
qualification, so that it is just as though the rules were defined originally
in `target_module`.

end::reference[] */
void clone_rules(std::tuple<std::string, std::string> source_target_modules);

/* tag::reference[]

== `b2::modules::call_in`

====
[horizontal]
Jam:: `rule call-in ( module-name ? : rule-name args * : * )`
{CPP}:: `list_ref call_in(value_ref module_name, std::tuple<value_ref, list_ref>
const lists & rest, bind::context_ref_ context_ref);`
====

Call the given rule locally in the given module. Use this for rules accepting
rule names as arguments, so that the passed rule may be invoked in the context
of the rule's caller (for example, if the rule accesses module globals or is a
local rule). Note that rules called this way may accept at most 18 parameters.

Example:

[source,jam]
----
rule filter ( f : values * )
{
	local m = [ CALLER_MODULE ] ;
	local result ;
	for v in $(values)
	{
		if [ modules.call-in $(m) : $(f) $(v) ]
		{
			result += $(v) ;
		}
	}
	return result ;
}
----

end::reference[] */
list_ref call_in(value_ref module_name,
	std::tuple<value_ref, list_ref> rule_name_a1,
	const lists & rest,
	bind::context_ref_ context_ref);

/* tag::reference[]

== `b2::modules::call_locally`

====
[horizontal]
Jam:: `rule call-locally ( qualified-rule-name args * : * )`
{CPP}:: `list_ref call_locally(std::tuple<value_ref, list_ref> rule_name_a1,
const lists & rest, const bind::context_ref_ & context_ref);`
====

Given a possibly qualified rule name and arguments, remove any initial module
qualification from the rule and invoke it in that module. If there is no
module qualification, the rule is invoked in the global module. Note that
rules called this way may accept at most 18 parameters.

end::reference[] */
list_ref call_locally(std::tuple<std::string, list_ref> rule_name_a1,
	const lists & rest,
	bind::context_ref_ context_ref);

/* tag::reference[]

== `b2::modules::run_tests`

====
[horizontal]
Jam:: `rule run-tests ( m )`
{CPP}:: `void run_tests(value_ref m, bind::context_ref_ context_ref);`
====

Runs internal B2 unit tests for the specified module. The module's
`__test__` rule is executed in its own module to eliminate any inadvertent
effects of testing module dependencies (such as assert) on the module itself.

end::reference[] */
void run_tests(value_ref m, bind::context_ref_ context_ref);

/* tag::reference[]

== `b2::modules::load`

====
[horizontal]
Jam:: `rule load ( module-name : filename ? : search * )`
{CPP}:: `void load(value_ref module_name, value_ref filename, list_cref
search, bind::context_ref_ context_ref);`
====

Load the indicated module if it is not already loaded.

`module-name`::
  Name of module to load.

`filename`::
  (partial) path to file; Defaults to `$(module-name).jam`

`search`::
  Directories in which to search for filename. Defaults to
  `$(BOOST_BUILD_PATH)`.

end::reference[] */
void load(value_ref module_name,
	value_ref filename,
	list_cref search,
	bind::context_ref_ context_ref);

/* tag::reference[]

== `b2::modules::import`

====
[horizontal]
Jam:: `rule import ( module-names + : rules-opt * : rename-opt * )`
{CPP}:: `void import(list_cref module_names, list_cref rules_opt, list_cref
rename_opt, bind::context_ref_ context_ref);`
====

Load the indicated module and import rule names into the current module. Any
members of `rules-opt` will be available without qualification in the caller's
module. Any members of `rename-opt` will be taken as the names of the rules in
the caller's module, in place of the names they have in the imported module.
If `rules-opt` = `*`, all rules from the indicated module are imported into the
caller's module. If `rename-opt` is supplied, it must have the same number of
elements as `rules-opt`.

NOTE: The `import` rule is available without qualification in all modules.

Examples:

[source,jam]
----
import path ;
import path : * ;
import path : join ;
import path : native make : native-path make-path ;
----

end::reference[] */
void import(list_cref module_names,
	list_cref rules_opt,
	list_cref rename_opt,
	bind::context_ref_ context_ref);

struct modules_module : b2::bind::module_<modules_module>
{
	const char * module_name = "modules";
	static const char * init_code;

	template <class Binder>
	void def(Binder & binder)
	{
		binder.def(&binding, "binding", "module_name" * _1)
			.def(&record_binding, "record-binding",
				"module_name" * _1 | "value" * _1)
			.def(&poke, "poke",
				"module_name" * _01 | "variables" * _1n | "value" * _n)
			.def(&peek, "peek", "module_name" * _01 | "variables" * _1n)
			.def(&clone_rules, "clone-rules",
				"source_module" * _1 + "target_module" * _1)
			.def(&call_in, "call-in",
				"module_name" * _01 | ("rule_name" * _1 + "a1" * _n)
					| "rest" * _r)
			.def(&call_locally, "call-locally",
				("rule_name" * _1 + "a1" * _n) | "rest" * _r)
			.def(&b2::paths::normalize_all, "normalize-paths", "paths" * _n)
			.def(&run_tests, "run-tests", "m" * _1)
			.def(&load, "load",
				"module-name" * _1 | "filename" * _01 | "search" * _n)
			.def(&import, "import",
				"module-names" * _1n | "rules-opt" * _n | "rename-opt" * _n);
		binder.eval(init_code);
		binder.loaded();
	}
};

}}} // namespace b2::jam::modules

#endif
