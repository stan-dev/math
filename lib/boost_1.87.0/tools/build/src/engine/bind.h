/*
Copyright 2019-2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_BIND_H
#define B2_BIND_H

#include "config.h"
#include <functional>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

/** tag::binder[]

= Binder

The B2 C++ native engine and system is reflected to various other languages
through bindings of classes and functions. This is accomplished through custom
C++ 11 reflection and non-intrusive declaration of a binding API.

end::binder[] */

namespace b2 { namespace bind {

/** tag::binder_type[]

== `b2::bind::type_`

Tag wrapper definition to specify a type to one of the binder specifications.
For example:

```
binder.def_class("system_info", type_<b2::system_info>());
```

end::binder_type[] */
template <class T>
struct type_
{
	typedef T type;
};

/** tag::binder_init[]

== `b2::bind::init_`

Tag wrapper definition to specify the types of arguments for init methods,
i.e. constructors, to a binder definition.
For example:

```
binder.def_class("system_info", type_<b2::system_info>())
	.def(init_<>());
```

end::binder_init[] */
template <class... Args>
struct init_
{};

// Forward declare..
template <int C>
struct arg_;
template <class... A>
struct args_;

/** tag::binder_param[]

A single parameter of an argument.

end::binder_param[] */
struct param_
{
	enum count_
	{
		one,
		any,
		many,
		optional,
        rest
	};

	// The symbolic name of this argument.
	const char * name = nullptr;

	// How many values this argument can accept.
	count_ count = one;

	param_() {}
	param_(const char * n, count_ c)
		: name(n)
		, count(c)
	{}
};

/** tag::binder_arg[]

A single argument's list of definitions.

end::binder_arg[] */
template <int C>
struct arg_
{
	enum
	{
		count = C
	};
	param_ args[C];
};

inline arg_<2> operator+(const param_ & a, const param_ & b)
{
	return { { a, b } };
}

inline arg_<3> operator+(const arg_<2> & a, const param_ & b)
{
	return { { a.args[0], a.args[1], b } };
}

inline arg_<4> operator+(const arg_<3> & a, const param_ & b)
{
	return { { a.args[0], a.args[1], a.args[2], b } };
}

template <int C>
arg_<C + 1> operator+(const arg_<C> & a, const param_ & b)
{
	arg_<C + 1> result;
	for (std::size_t i = 0; i < C - 1; ++i) result.args[i] = a.args[i];
	result.args[C] = b;
	return result;
}

/** tag::binder_args[]
end::binder_args[] */
template <class... A>
struct args_
{
	enum
	{
		count = sizeof...(A)
	};
	std::tuple<A...> arg;
};

template <class... A, int BC>
auto operator|(const args_<A...> & a, const arg_<BC> & b)
	-> args_<A..., arg_<BC>>
{
	return { { std::tuple_cat(a.arg, std::make_tuple(b)) } };
}

template <int AC, int BC>
auto operator|(const arg_<AC> & a, const arg_<BC> & b)
	-> args_<arg_<AC>, arg_<BC>>
{
	return { { std::make_tuple(a, b) } };
}

inline auto operator|(const param_ & a, const param_ & b)
	-> args_<arg_<1>, arg_<1>>
{
	return arg_<1> { { a } } | arg_<1> { { b } };
}

template <int AC>
auto operator|(const arg_<AC> & a, const param_ & b) -> args_<arg_<AC>, arg_<1>>
{
	return a | arg_<1> { { b } };
}

template <class... A>
auto operator|(const args_<A...> & a, const param_ & b) -> args_<A..., arg_<1>>
{
	return a | arg_<1> { { b } };
}

template <int BC>
auto operator|(const param_ & a, const arg_<BC> & b) -> args_<arg_<1>, arg_<BC>>
{
	return arg_<1> { { a } } | b;
}

/** tag::binder_module[]

== `b2::bind::module_`

The base type for a module binding. Modules are a collection of class and other
binding declarations. A module is specified as a concrete class with at minimum
a `module_name` and a definition method
(`template <class Binder> void def(Binder & binder)`).
For example:

```
struct sysinfo_module
	: b2::bind::module_<sysinfo_module>
{
	const char *module_name = "sysinfo";

	template <class Binder>
	void def(Binder & binder)
	{
		// ...
	}
};
```

end::binder_module[] */
template <class Module>
struct module_
{
	// Alias shorthand for `b2::bind::type_` to avoid namespace qualification.
	template <class T>
	using type_ = ::b2::bind::type_<T>;

	// Alias shorthand for `b2::bind::init_` to avoid namespace qualification.
	template <class... A>
	using init_ = ::b2::bind::init_<A...>;

	// Get the sub-typed, from CRTP, `*this` reference.
	Module & self() { return *static_cast<Module *>(this); }
};

/** tag::binder_class[]

== `b2::bind::class_`

The `class_` type provides for defining class level members, i.e. methods. The
`class_` type is returned by the `binder_::def_class(..)` to give a quick, and
scoped, manner to define complete class bindings.
For example:

```
binder.def_class("system_info", type_<b2::system_info>())
	.def(init_<>())
	.def("cpu_core_count", &b2::system_info::cpu_core_count)
	.def("cpu_thread_count", &b2::system_info::cpu_thread_count);
```

end::binder_class[] */
template <class Class, class Binder>
struct class_
{
	Binder & binder;

	class_(const char * name, Binder & binder)
		: binder(binder)
	{
		name_ref() = name;
	}

	static const char * name() { return name_ref().c_str(); }

	/** tag::binder_class[]

	=== `b2::bind::class_::def(init_<Args...> init_args, ...) -> class_ &`

	Defines an init method, i.e. constructor, with the given `init_args` types.
	The given `init_args` will match an appropriate constructor. And the
	constructor will be called from the bound language as appropriate for the
	language. For example for `jam` it will call the `__init__` rule.

	end::binder_class[] */
	template <class... Args>
	class_ & def(init_<Args...> init_args)
	{
		// Forward to the language specific binder.
		binder.def_init(this->name(), (Class *)nullptr, init_args, args_<> {});
		return *this;
	}
	template <class... Args, class... A>
	class_ & def(init_<Args...> init_args, args_<A...> args)
	{
		// Forward to the language specific binder.
		binder.def_init(this->name(), (Class *)nullptr, init_args, args);
		return *this;
	}

	/** tag::binder_class[]

	=== `b2::bind::class_::def(F function, const char *name, args_<A...> args)
	-> class_ &`

	Defines a method which is bound with the given `name` which calls the given
	`function`.

	end::binder_class[] */
	template <class F, class... A>
	class_ & def(F function, const char * name, args_<A...> args)
	{
		// Forward to the language specific binder.
		binder.def_method(this->name(), name, args, function);
		return *this;
	}
	template <class F, int C>
	class_ & def(F function, const char * name, arg_<C> args)
	{
		return this->def(
			function, name, args_<arg_<C>> { { std::make_tuple(args) } });
	}
	template <class F>
	class_ & def(F function, const char * name, param_ args)
	{
		return this->def(function, name,
			args_<arg_<1>> { { std::make_tuple(arg_<1> { { args } }) } });
	}
	template <class F>
	class_ & def(F function, const char * name)
	{
		return this->def(function, name, args_<> {});
	}

	private:
	static std::string & name_ref()
	{
		static std::string name;
		return name;
	}
};

/* tag::binder[]

Interface that converts between binder values and C++ values.

end::binder[] */
template <class Binder, class CxxValue, class BindValue>
struct converter_
{
	static BindValue to_bind_value(const CxxValue &);
	static CxxValue from_bind_value(BindValue);
};

struct context_
{};

struct context_ref_
{
	explicit context_ref_(context_ & c)
		: context(&c)
	{}
	context_ref_() = default;
	context_ref_(const context_ref_ &) = default;
	context_ref_(context_ref_ && o)
	{
		context = o.context;
		o.context = nullptr;
	}
	context_ref_ & operator=(const context_ref_ &) = default;
	context_ref_ & operator=(context_ & c)
	{
		context = &c;
		return *this;
	}

	template <typename C>
	C & get() const
	{
		return static_cast<C &>(*context);
	}

	private:
	context_ * context = nullptr;
};

struct binder_interface_
{
	context_ref_ context_ref;
};

/** tag::binder_binder[]

== `b2::bind::binder_`

Interface for a language agnostic binder. This is passed to the module
definition function (`template <class B> def(B & binder)`). What is passed
is a language specific subclass that will generate the needed bindings for
that language.

end::binder_binder[] */
template <class Binder>
struct binder_ : binder_interface_
{
	/** tag::binder_binder[]

	=== `b2::bind::binder_::def_class`

	Declares the definition of a class, given in the `type_` wrapper, in the
	module. Optional base class types can be given as additional arguments.
	The base classes need to be bound before this subclass is bound.

	end::binder_binder[] */
	template <class Class, class... Bases>
	class_<Class, Binder> def_class(
		const char * name, type_<Class> class_type, type_<Bases>... bases)
	{
		class_<Class, Binder> class_def { name, self() };
		self().bind_class(current_module_name, name, class_type, bases...);
		return class_def;
	}

	/** tag::binder_binder[]

	=== `b2::bind::binder_::def_class`

	=== `b2::bind::binder_::def(F function, const char *name, ...) -> binder_ &`

	Defines a function which is bound with the given `name` which calls the
	given `function`.

	end::binder_binder[] */
	template <class F, class... A>
	binder_ & def(F function, const char * name, args_<A...> args)
	{
		// Forward to the language specific binder.
		self().def_function(name, args, function);
		return *this;
	}
	template <class F, int C>
	binder_ & def(F function, const char * name, arg_<C> args)
	{
		return this->def(
			function, name, args_<arg_<C>> { { std::make_tuple(args) } });
	}
	template <class F>
	binder_ & def(F function, const char * name, param_ args)
	{
		return this->def(function, name,
			args_<arg_<1>> { { std::make_tuple(arg_<1> { { args } }) } });
	}
	template <class F>
	binder_ & def(F function, const char * name)
	{
		return this->def(function, name, args_<> {});
	}

	binder_ & eval(const char * data)
	{
		self().eval_data(current_module_name, data);
		return *this;
	}

	void loaded() { self().set_loaded(current_module_name); }

	// Internal..

	// Returns the subclass reference to this binder.
	Binder & self() { return *static_cast<Binder *>(this); }

	// Binds the given native module declarations. This calls the subclass'
	// `bind_module(module_name)` method to do any binding work for the module.
	// And then calls the `def(binder)` method on the module to run through the
	// definitions of the module for this binder.
	template <class Module>
	Binder & bind(Module m)
	{
		current_module_name = m.module_name;
		self().bind_module(current_module_name);
		m.def(self());
		return self();
	}

	// Respond to a method definition of a class. This calls the subclass
	// method `bind_method(module_name, class_name, method_name, function)`.
	template <class Function, class... A>
	void def_method(const char * class_name,
		const char * method_name,
		args_<A...> args,
		Function f)
	{
		self().bind_method(
			current_module_name, class_name, method_name, args, f);
	}

	// Respond to a constructor definition of a class. This calls the subclass
	// method `bind_init(module_name, class_name, class_nullptr, init)`.
	template <class Class, class Init, class... A>
	void def_init(const char * class_name, Class * c, Init i, args_<A...> args)
	{
		self().bind_init(current_module_name, class_name, c, i, args);
	}

	// Respond to a function definition of a module. This calls the subclass
	// method `bind_function(module_name, function_name, function)`.
	template <class Function, class... A>
	void def_function(const char * function_name, args_<A...> args, Function f)
	{
		self().bind_function(current_module_name, function_name, args, f);
	}

	// Generic, shim, to convert from a C++ value to a binding specific value.
	// Forwards to the `converter_` template specialization.
	template <class BindValue, class CxxValue>
	static BindValue convert_to_bind_value(const CxxValue & source)
	{
		return converter_<Binder, CxxValue, BindValue>::to_bind_value(source);
	}

	// Generic, shim, to convert from a binding specific value to a C++ value.
	// Forwards to the `converter_` template specialization.
	template <class CxxValue, class BindValue>
	static CxxValue convert_from_bind_value(BindValue source)
	// static CxxValue convert_from_bind_value(const BindValue &source)
	{
		return converter_<Binder, CxxValue, BindValue>::from_bind_value(source);
	}

	protected:
	const char * current_module_name = nullptr;
};

}} // namespace b2::bind

namespace b2 {

enum bind_param_count_one : bool
{
	_1 = true
};
enum bind_param_count_any : bool
{
	_n = true
};
enum bind_param_count_many : bool
{
	_1n = true
};
enum bind_param_count_optional : bool
{
	_01 = true
};
enum bind_param_count_rest : bool
{
	_r = true
};

inline bind::param_ operator*(const char * name, bind_param_count_one)
{
	return bind::param_ { name, bind::param_::one };
}

inline bind::param_ operator*(const char * name, bind_param_count_any)
{
	return bind::param_ { name, bind::param_::any };
}

inline bind::param_ operator*(const char * name, bind_param_count_many)
{
	return bind::param_ { name, bind::param_::many };
}

inline bind::param_ operator*(const char * name, bind_param_count_optional)
{
	return bind::param_ { name, bind::param_::optional };
}

inline bind::param_ operator*(const char * name, bind_param_count_rest)
{
	return bind::param_ { name, bind::param_::rest };
}

} // namespace b2

#endif
