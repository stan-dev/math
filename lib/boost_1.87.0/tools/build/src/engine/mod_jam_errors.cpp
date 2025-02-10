/*
Copyright 2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "mod_jam_errors.h"

#include "bindjam.h"
#include "startup.h"
#include "variable.h"

namespace b2 { namespace jam { namespace errors {

namespace {

inline lists & last_error()
{
	static lists e;
	return e;
}

inline std::size_t & disabled()
{
	static std::size_t d = 0;
	return d;
}

inline bool is_disabled() { return disabled() > 0; }
inline void set_disabled() { disabled() += 1; }
inline bool clear_disabled() { return (disabled() -= 1) > 0; }

} // namespace

void backtrace(int skip_frames,
	const std::string & prefix,
	const lists & messages,
	bind::context_ref_ context_ref)
{
	bool user_modules_only
		= !variable("errors", ".user-modules-only").get().empty();

	FRAME * frame = context_ref.get<jam_context>().frame;
	bool print_args = true;
	do
	{
		if (!user_modules_only && skip_frames-- > 0) continue;
		if (user_modules_only)
			frame = frame->module->user_module ? frame : frame->prev_user;

		std::string info;
		if (frame)
		{
			if (user_modules_only) info += prefix + " at ";
			info += frame->file ? object_str(frame->file) : "(builtin)";
			info += ":";
			info += std::to_string(frame->line);
			info += ": in ";
			info += frame->rulename;
			if (frame->module->name)
			{
				info += " from module ";
				info += object_str(frame->module->name);
			}
			info += "\n";
		}
		if (!info.empty()) out_puts(info.c_str());

		// The first time through, print each argument on a separate line.
		if (print_args)
		{
			print_args = false;
			auto print_arg = [&prefix](const list_cref & arg) {
				if (arg.empty()) return;
				if (!prefix.empty()) out_printf("%s ", prefix.c_str());
				list_print(*arg);
				out_putc('\n');
			};
			for (lists::size_type i = 0; i < messages.length(); ++i)
				print_arg(messages[i]);
		}
	}
	while (!user_modules_only && ((frame = frame->prev) != nullptr));
}

void backtrace(std::tuple<int, std::string, list_ref> skip_prefix_messages,
	const lists & rest,
	bind::context_ref_ context_ref)
{
	lists messages;
	messages.push_back(std::move(std::get<2>(skip_prefix_messages)));
	messages.append(rest);
	backtrace(std::get<0>(skip_prefix_messages),
		std::get<1>(skip_prefix_messages), messages, context_ref);
}

void error_skip_frames(
	int skip, const lists & messages, bind::context_ref_ context_ref)
{
	if (!is_disabled())
	{
		backtrace(skip, "error:", messages, context_ref);
		clean_exit(exit_result::failure);
	}
	else if (last_error().empty())
	{
		last_error().append(messages);
	}
}

void error_skip_frames(std::tuple<int, list_ref> skip_messages,
	const lists & rest,
	bind::context_ref_ context_ref)
{
	lists messages;
	messages.push_back(std::move(std::get<1>(skip_messages)));
	messages.append(rest);
	error_skip_frames(std::get<0>(skip_messages), messages, context_ref);
}

void error_try()
{
	set_disabled();
	last_error().clear();
}

void error_catch(const lists & messages, bind::context_ref_ context_ref)
{
	clear_disabled();
	if (last_error().empty())
		error_skip_frames(1,
			lists() | list_ref("expected an error, but none occurred"),
			context_ref);
	else
	{
		for (lists::size_type i = 0; i < messages.length(); ++i)
		{
			if (i < last_error().length())
			{
				std::string expected;
				for (auto v : messages[i])
				{
					if (!expected.empty()) expected += " ";
					expected += v->str();
				}
				std::string got;
				for (auto v : last_error()[i])
				{
					if (!got.empty()) got += " ";
					got += v->str();
				}
				if (expected != got)
				{
					last_error().clear();
					error_skip_frames(1,
						lists()
							| list_ref("expected \"" + expected
								+ "\" in argument " + std::to_string(i + 1)
								+ " of error")
							| list_ref("got \"" + got + "\" instead"),
						context_ref);
				}
			}
		}
	}
}

void error(const lists & rest, bind::context_ref_ context_ref)
{
	if (!variable("errors", ".no-error-backtrace").get().empty())
	{
		// Print each argument on a separate line.
		for (lists::size_type i = 0; i < rest.length(); ++i)
		{
			if (i == 0) out_puts("error: ");
			list_print(*rest[i]);
			out_putc('\n');
		}
		clean_exit(exit_result::failure);
	}
	else
		error_skip_frames(1, rest, context_ref);
}

void user_error(const lists & rest, bind::context_ref_ context_ref)
{
	variable("errors", ".user-modules-only") = list_ref() + "1";
	error_skip_frames(1, rest, context_ref);
}

void warning(const lists & rest, bind::context_ref_ context_ref)
{
	backtrace(1, "warning:", rest, context_ref);
}

list_ref lol_to_list(const lists & rest)
{
	list_ref result;
	lists::size_type size = rest.length();
	while (size > 0 && rest[size - 1].empty()) --size;
	for (lists::size_type i = 0; i < size; ++i)
	{
		if (i != 0) result.push_back(":");
		for (auto l : rest[i])
		{
			result.push_back("\"" + std::string(l->str()) + "\"");
		}
	}
	return result;
}

list_ref nearest_user_location(bind::context_ref_ context_ref)
{
	FRAME * frame = context_ref.get<jam_context>().frame;
	frame = frame->module->user_module ? frame : frame->prev_user;
	if (!frame) return list_ref();
	return list_ref(
		std::string(frame->file ? object_str(frame->file) : "(builtin)") + ":"
		+ std::to_string(frame->line));
}
/*
Jam init code is:

# Copyright 2003 Dave Abrahams
# Copyright 2004 Vladimir Prus
*/
const char * errors_module::init_code = R"jam(

if --no-error-backtrace in [ modules.peek : ARGV ]
{
    .no-error-backtrace = true ;
}

rule __test__ ( )
{
    import assert ;

    assert.result \"a\" \"b\" \"c\" ":" \"d\" : lol->list a b c : d ; 

    # Show that we can correctly catch an expected error.
    try ;
    {
        error an error occurred : somewhere ;
    }
    catch an error occurred : somewhere ;

    # Show that unexpected errors generate real errors.
    try ;
    {
        try ;
        {
            error an error occurred : somewhere ;
        }
        catch an error occurred : nowhere ;
    }
    catch expected \"nowhere\" in argument 2 of error ;

    # Show that not catching an error where one was expected is an error.
    try ;
    {
        try ;
        {
        }
        catch ;
    }
    catch expected an error, but none occurred ;
}

)jam";

}}} // namespace b2::jam::errors
