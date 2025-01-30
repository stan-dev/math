/*
Copyright 2019-2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_BINDJAM_H
#define B2_BINDJAM_H

#include "bind.h"
#include "frames.h"
#include "value.h"

namespace b2 { namespace jam {

/*
Context for Jam calls.
*/
struct jam_context : bind::context_
{
	FRAME * frame = nullptr;

	jam_context(FRAME * f)
		: frame(f)
	{}

	inline bind::context_ref_ ref() { return bind::context_ref_(*this); }
};

/*
Jam language C++ interface binder.
*/
struct jam_binder : b2::bind::binder_<jam_binder>
{
	/*
	Bind the named module to Jam. All this does for Jam is define the module
	scope in Jam.
	*/
	void bind_module(const char * module_name);

	/*
	Register the class "meta-class". Which in the case of jam is a module
	with the rules and variables that gets imported and copied into the
	specific instance module as needed.
	*/
	template <class Class, class... Bases>
	void bind_class(const char * module_name,
		const char * class_name,
		::b2::bind::type_<Class>,
		::b2::bind::type_<Bases>...);

	/*
	Bind the init ctor function `i` as a constructor for the class.
	*/
	template <class Class, class Init, class... A>
	void bind_init(const char * module_name,
		const char * class_name,
		Class * c,
		Init i,
		::b2::bind::args_<A...> args);

	/*
	Bind the given method of a class. THe function can be any invocable whose
	interface will be reflected to Jam.
	*/
	template <class Function, class... A>
	void bind_method(const char * module_name,
		const char * class_name,
		const char * method_name,
		::b2::bind::args_<A...> args,
		Function f);

	/*
	Bind a module scoped plain function, i.e. Jam rule.
	*/
	template <class Function, class... A>
	void bind_function(const char * module_name,
		const char * function_name,
		::b2::bind::args_<A...> args,
		Function f);

	/*
	Evaluate the Jam program in `data` within the `modules_name` scope.
	*/
	void eval_data(const char * module_name, const char * data);

	/*
	Mark the module as loaded. Future attempts at loading, with `module.load`,
	will fail (as noop).
	*/
	void set_loaded(const char * module_name);
};

/*
Marshaling of arguments and results end up calling one of `from_jam` or
`to_jam` functions. These are specialized on some core types as needed:
`std::string`, `unsigned int`, `int`, `bool`
*/

template <class Value>
Value from_jam(value_ptr jam_value);

template <class Value>
value_ptr to_jam(Value value);

/*
Binds all the declared C++ interfaces to Jam equivalents.
*/
template <typename F>
void bind_jam(F * frame);

}} // namespace b2::jam

#endif
