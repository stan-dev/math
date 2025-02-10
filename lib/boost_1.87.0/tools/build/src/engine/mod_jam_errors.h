/*
Copyright 2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_MOD_JAM_ERRORS_H
#define B2_MOD_JAM_ERRORS_H

#include "config.h"

#include "bind.h"
#include "lists.h"
#include "value.h"

#include <tuple>

/* tag::reference[]

[[b2.reference.modules.errors]]
= `errors` module.

end::reference[] */

namespace b2 { namespace jam { namespace errors {

/* tag::reference[]

== `b2::jam::errors::backtrace`

====
[horizontal]
Jam:: `rule backtrace ( skip-frames prefix messages * : * )`
{CPP}:: `void backtrace(std::tuple<int, std::string, list_ref>
skip_prefix_messages, const lists & rest, bind::context_ref_ context_ref);`
====

Print a stack backtrace leading to this rule's caller. Each argument
represents a line of output to be printed after the first line of the
backtrace.

end::reference[] */
void backtrace(std::tuple<int, std::string, list_ref> skip_prefix_messages,
	const lists & rest,
	bind::context_ref_ context_ref);

/* tag::reference[]

== `b2::jam::errors::error_skip_frames`

====
[horizontal]
Jam:: `rule error-skip-frames ( skip-frames messages * : * )`
{CPP}:: `void error_skip_frames(std::tuple<int, list_ref> skip_messages, const
lists & rest, bind::context_ref_ context_ref);`
====

end::reference[] */
void error_skip_frames(std::tuple<int, list_ref> skip_messages,
	const lists & rest,
	bind::context_ref_ context_ref);

/* tag::reference[]

== `b2::jam::errors::try-catch`

====
[horizontal]
Jam:: `rule error-skip-frames ( skip-frames messages * : * )`
{CPP}:: `void error_skip_frames(std::tuple<int, list_ref> skip_messages, const
lists & rest, bind::context_ref_ context_ref);`
====

This is not really an exception-handling mechanism, but it does allow us to
perform some error-checking on our error-checking. Errors are suppressed
after a try, and the first one is recorded. Use catch to check that the
error message matched expectations.

====
[horizontal]
Jam:: `rule try ( )`
{CPP}:: `void error_try();`
====

Begin looking for error messages.

====
[horizontal]
Jam:: `rule catch ( messages * : * )`
{CPP}:: `void error_catch(const lists & rest, bind::context_ref_ context_ref);`
====

Stop looking for error messages; generate an error if an argument of
messages is not found in the corresponding argument in the error call.

end::reference[] */
void error_try();
void error_catch(const lists & rest, bind::context_ref_ context_ref);

/* tag::reference[]

== `b2::jam::errors::error`

====
[horizontal]
Jam:: `rule error ( messages * : * )`
{CPP}:: `void error(const lists & rest, bind::context_ref_ context_ref);`
====

Print an error message with a stack backtrace and exit.

end::reference[] */
void error(const lists & rest, bind::context_ref_ context_ref);

/* tag::reference[]

== `b2::jam::errors::user_error`

====
[horizontal]
Jam:: `rule user-error ( messages * : * )`
{CPP}:: `void user_error(const lists & rest, bind::context_ref_ context_ref);`
====

Same as 'error', but the generated backtrace will include only user files.
end::reference[] */
void user_error(const lists & rest, bind::context_ref_ context_ref);

/* tag::reference[]

== `b2::jam::errors::warning`

====
[horizontal]
Jam:: `rule warning ( messages * : * )`
{CPP}:: `void warning(const lists & rest, bind::context_ref_ context_ref);`
====

Print a warning message with a stack backtrace and exit.

end::reference[] */
void warning(const lists & rest, bind::context_ref_ context_ref);

/* tag::reference[]

== `b2::jam::errors::lol_to_list`

====
[horizontal]
Jam:: `rule lol->list ( * )`
{CPP}:: `list_ref lol_to_list(const lists & rest);`
====

Convert an arbitrary argument list into a list with ":" separators and quoted
elements representing the same information. This is mostly useful for
formatting descriptions of arguments with which a rule was called when
reporting an error.

end::reference[] */
list_ref lol_to_list(const lists & rest);

/* tag::reference[]

== `b2::jam::errors::nearest_user_location`

====
[horizontal]
Jam:: `rule nearest-user-location ( )`
{CPP}:: `list_ref nearest_user_location(bind::context_ref_ context_ref);`
====

Return the file:line for the nearest entry in backtrace which correspond to a
user module.

end::reference[] */
list_ref nearest_user_location(bind::context_ref_ context_ref);

struct errors_module : b2::bind::module_<errors_module>
{
	const char * module_name = "errors";
	static const char * init_code;

	template <class Binder>
	void def(Binder & binder)
	{
		binder
			.def(&backtrace, "backtrace",
				("skip-frames" * _1 + "prefix" * _1 + "messages" * _n)
					| ("rest" * _r))
			.def(&error_skip_frames, "error-skip-frames",
				("skip-frames" * _1 + "messages" * _n) | ("rest" * _r))
			.def(&error_try, "try")
			.def(&error_catch, "catch", "messages" * _r)
			.def(&error, "error", "messages" * _r)
			.def(&user_error, "user-error", "messages" * _r)
			.def(&warning, "warning", "messages" * _r)
			.def(&nearest_user_location, "nearest-user-location")
			.def(&lol_to_list, "lol->list", "rest" * _r);
		binder.eval(init_code);
		binder.loaded();
	}
};

}}} // namespace b2::jam::errors

#endif
