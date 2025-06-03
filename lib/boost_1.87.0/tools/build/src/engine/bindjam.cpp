/*
Copyright 2019-2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "bindjam.h"

#include "class.h"
#include "function.h"
#include "hash.h"
#include "lists.h"
#include "modules.h"
#include "mp.h"
#include "native.h"
#include "optval.h"
#include "output.h"
#include "parse.h"
#include "rules.h"
#include "types.h"
#include "value.h"
#include "variable.h"

#include "mod_command_db.h"
#include "mod_db.h"
#include "mod_jam_builtin.h"
#include "mod_jam_class.h"
#include "mod_jam_errors.h"
#include "mod_jam_modules.h"
#include "mod_path.h"
#include "mod_regex.h"
#include "mod_set.h"
#include "mod_string.h"
#include "mod_sysinfo.h"
#include "mod_version.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace b2 { namespace jam {

namespace {
inline const char * mod_cstr(const std::string & n)
{
	return n.empty() ? nullptr : n.c_str();
}
inline std::string mod_str(const char * n) { return (n == nullptr) ? "" : n; }
inline module_ptr mod(const char * name)
{
	if (name == nullptr || name[0] == '\0') return root_module();
	return bindmodule(value_ref(name));
}
inline module_ptr mod(const std::string & name) { return mod(mod_cstr(name)); }
} // namespace
/*
Basic core types to marshal..
*/

template <>
string_t from_jam<string_t>(value_ptr jv)
{
	return object_str(jv);
}
template <>
value_ptr to_jam(string_t v)
{
	return object_new(v.c_str());
}

template <>
uint_t from_jam<uint_t>(value_ptr jv)
{
	return std::stoul(from_jam<string_t>(jv));
}
template <>
value_ptr to_jam(uint_t v)
{
	return to_jam(std::to_string(v));
}

template <>
int_t from_jam<int_t>(value_ptr jv)
{
	return std::stoi(from_jam<string_t>(jv));
}
template <>
value_ptr to_jam(int_t v)
{
	return to_jam(std::to_string(v));
}

template <>
bool from_jam<bool>(value_ptr jv)
{
	return object_str(jv) != nullptr && object_str(jv)[0] != '\0';
}
template <>
value_ptr to_jam(bool v)
{
	return object_new(v ? "1" : "");
}

template <>
value_ref from_jam<value_ref>(value_ptr v)
{
	return value_ref(v);
}

// General marshaling of one jam value list. Default converts the first item
// to/from the list.
template <class CxxValue>
struct converter
{
	static list_ref to_bind_value(const CxxValue & cxx_value)
	{
		return list_ref(to_jam(cxx_value));
	}
	static CxxValue from_bind_value(
		list_cref::iterator & i, list_cref::iterator e)
	{
		return (i == e) ? CxxValue() : from_jam<CxxValue>(*(i++));
	}
};

// Marshaling specialization for vector container.
// TODO: Generalize to more containers.
template <class Value>
struct converter<std::vector<Value>>
{
	using CxxValue = std::vector<Value>;
	static list_ref to_bind_value(const CxxValue & cxx_value)
	{
		list_ref result;
		for (auto & v : cxx_value) result.push_back(to_jam(v));
		return result;
	}
	static CxxValue from_bind_value(
		list_cref::iterator & i, list_cref::iterator e)
	{
		CxxValue result;
		for (; i != e; ++i) result.emplace_back(from_jam<Value>(*i));
		return result;
	}
};

// Marshaling for optval container.
template <class ValueType>
struct converter<optval<ValueType>>
{
	using CxxValue = optval<ValueType>;
	static list_ref to_bind_value(const CxxValue & cxx_value)
	{
		list_ref result;
		if (cxx_value.has_value())
			result.push_back(to_jam(static_cast<ValueType>(cxx_value)));
		return result;
	}
	static CxxValue from_bind_value(
		list_cref::iterator & i, list_cref::iterator e)
	{
		return (i == e) ? CxxValue() : from_jam<ValueType>(*(i++));
	}
};

//
template <>
struct converter<list_ref>
{
	static list_ref to_bind_value(const list_ref & cxx_value)
	{
		return cxx_value;
	}
	static list_ref from_bind_value(
		list_cref::iterator & i, list_cref::iterator e)
	{
		list_ref result(i, e);
		i = e;
		return result;
	}
};

// Marshaling for value_ref container.
template <>
struct converter<value_ref>
{
	static list_ref to_bind_value(const value_ref & cxx_value)
	{
		list_ref result;
		if (cxx_value.has_value()) result.push_back(cxx_value);
		return result;
	}
	static value_ref from_bind_value(
		list_cref::iterator & i, list_cref::iterator e)
	{
		return (i == e) ? value_ref() : value_ref(*(i++));
	}
};

// Marshaling of tuples as multi-param arguments.
template <class... ValueTypes>
struct converter<std::tuple<ValueTypes...>>
{
	using CxxValue = std::tuple<ValueTypes...>;

	static list_ref to_bind_value(const CxxValue & cxx_value)
	{
		return to_bind_value(cxx_value,
			mp::make_index_sequence<std::tuple_size<CxxValue>::value> {});
	}
	template <std::size_t... I>
	static list_ref to_bind_value(
		const CxxValue & cxx_value, mp::index_sequence<I...>)
	{
		using std::get;
		value_ptr jam_val[] { nullptr, to_jam(get<I>(cxx_value))... };
		list_ref result;
		for (int i = 1; i <= std::tuple_size<CxxValue>::value; ++i)
			result.push_back(jam_val[i]);
		return result;
	}

	static CxxValue from_bind_value(
		list_cref::iterator & i, list_cref::iterator e)
	{
		return from_bind_value(
			i, e, mp::make_index_sequence<std::tuple_size<CxxValue>::value> {});
	}
	template <std::size_t... I>
	static CxxValue from_bind_value(list_cref::iterator & i,
		list_cref::iterator e,
		mp::index_sequence<I...>)
	{
		CxxValue result;
		int _[] { 0,
			((void)(std::get<I>(result)
				 = converter<typename std::tuple_element<I,
					 CxxValue>::type>::from_bind_value(i, e)),
				0)... };
		(void)_;
		return result;
	}
};

}} // namespace b2::jam

namespace b2 { namespace bind {

using namespace b2::jam;

// Direct, marshalling for list_cref container.
template <>
struct converter_<b2::jam::jam_binder, list_cref, LIST *>
{
	static LIST * to_bind_value(const list_cref & cpp_value)
	{
		return *cpp_value;
	}
	static list_cref from_bind_value(LIST * jam_value)
	{
		return list_cref(jam_value);
	}
};

// Direct, marshalling for list_ref container.
template <>
struct converter_<b2::jam::jam_binder, list_ref, LIST *>
{
	static LIST * to_bind_value(const list_ref & cpp_value)
	{
		return list_ref(cpp_value).release();
	}
	static list_ref from_bind_value(LIST * jam_value)
	{
		// TODO: Can we optimize this to not do a list copy?
		return list_ref(jam_value);
	}
};

template <class CxxValue>
struct converter_<b2::jam::jam_binder, CxxValue, LIST *>
{
	static LIST * to_bind_value(const CxxValue & cxx_value)
	{
		return b2::jam::converter<CxxValue>::to_bind_value(cxx_value).release();
	}
	static CxxValue from_bind_value(LIST * jam_value)
	{
		list_cref v(jam_value);
		auto i = v.begin();
		auto e = v.end();
		return b2::jam::converter<CxxValue>::from_bind_value(i, e);
	}
};

}} // namespace b2::bind

namespace b2 { namespace jam {

struct jam_cxx_self
{
	// Get the hidden C++ self instance for the jam class instance.
	template <class Class>
	static Class * get(FRAME * frame)
	{
		value_ref cxx_self_name(value::make("__cxx_self__"));
		value_ref cxx_self_value(
			list_front(var_get(frame->module, cxx_self_name)));
		Class * self = static_cast<Class *>(cxx_self_value->as_object());
		return self;
	}

	// Set the hidden C++ self instance for the jam class instance.
	template <class Class>
	static void set(FRAME * frame, Class * self)
	{
		value_ref cxx_self_name(value::make("__cxx_self__"));
		value_ref cxx_self_value(value::make(self));
		var_set(
			frame->module, cxx_self_name, list_new(cxx_self_value), VAR_SET);
	}
};

// Marshals arguments from Jam function FRAME to C++ tuple..

template <std::size_t N, class ArgType, class Enable = void>
struct jam_marshal_arg
{
	template <class ArgsTuple>
	static void convert(ArgsTuple & to_args, jam_context & context)
	{
		LIST * arg_value = lol_get(context.frame->args, N);
		if (arg_value != L0)
			std::get<N>(to_args)
				= jam_binder::convert_from_bind_value<ArgType>(arg_value);
	}
};

template <std::size_t N>
struct jam_marshal_arg<N, bind::context_ref_>
{
	template <class ArgsTuple>
	static void convert(ArgsTuple & to_args, jam_context & context)
	{
		std::get<N>(to_args) = bind::context_ref_(context);
	}
};

template <std::size_t N>
struct jam_marshal_arg<N, lists>
{
	template <class ArgsTuple>
	static void convert(ArgsTuple & to_args, jam_context & context)
	{
		lists & to_arg = std::get<N>(to_args);
		for (int32_t i = N; i < context.frame->args->count; ++i)
			to_arg | lol_get(context.frame->args, i);
	}
};

template <class ArgsTuple, std::size_t N, class Enable = void>
struct jam_marshal_args
{
	static void convert(ArgsTuple & to_args, jam_context & context) {}
};
template <class ArgsTuple, std::size_t N>
struct jam_marshal_args<ArgsTuple,
	N,
	typename std::enable_if<(N < std::tuple_size<ArgsTuple>::value)>::type>
{
	static void convert(ArgsTuple & to_args, jam_context & context)
	{
		using arg_type = typename std::tuple_element<N, ArgsTuple>::type;
		jam_marshal_arg<N, arg_type>::convert(to_args, context);

		if (N < std::tuple_size<ArgsTuple>::value)
		{
			jam_marshal_args<ArgsTuple, N + 1>::convert(to_args, context);
		}
	}
};

// Bound init/constructor function that forwards from jam __init__ to C++.
template <class Call, class Class, class... Args>
static LIST * jam_call_init(
	FRAME * frame, int flags, Call cxx_call, Class * (*)(Args...))
{
	using namespace mp;
	try
	{
		// Marshal arguments from frame->args.
		using ArgsTuple = std::tuple<typename remove_cvref<Args>::type...>;
		ArgsTuple args;
		jam_context context { frame };
		jam_marshal_args<ArgsTuple, 0>::convert(args, context);
		// Construct the class instance, and set the hidden jam instance var
		// to keep track of it.
		Class * self = invoke<Class *>(cxx_call, args,
			make_index_sequence<std::tuple_size<ArgsTuple>::value> {});
		jam_cxx_self::set(frame, self);
		// Nothing to return from an init.
		return L0;
	}
	catch (const std::exception &)
	{
		return L0;
	}
}

// Bound plain function that forwards from jam to C++ with a return
// converted back to jam.
template <class Call,
	class Class,
	class... Args,
	class Return,
	typename std::enable_if<!std::is_void<Return>::value, int>::type = 0>
static LIST * jam_call_method(
	FRAME * frame, int flags, Call cxx_call, Return (*)(Class *, Args...))
{
	using namespace mp;
	try
	{
		// Marshal arguments from frame->args.
		using ArgsTuple = std::tuple<typename remove_cvref<Args>::type...>;
		ArgsTuple args;
		jam_context context { frame };
		jam_marshal_args<ArgsTuple, 0>::convert(args, context);
		// Invoke call, and return marshaled result.
		return jam_binder::convert_to_bind_value<LIST *>(
			invoke<Return>(cxx_call,
				std::tuple_cat(
					std::make_tuple(jam_cxx_self::get<Class>(frame)), args),
				make_index_sequence<1 + std::tuple_size<ArgsTuple>::value> {}));
	}
	catch (const std::exception &)
	{
		return L0;
	}
}

// Bound plain function that forwards from jam to C++ with a void return.
template <class Call,
	class Class,
	class... Args,
	class Return,
	typename std::enable_if<std::is_void<Return>::value, int>::type = 0>
static LIST * jam_call_method(
	FRAME * frame, int flags, Call cxx_call, Return (*)(Class *, Args...))
{
	using namespace mp;
	try
	{
		// Marshal arguments from frame->args.
		using ArgsTuple = std::tuple<typename remove_cvref<Args>::type...>;
		ArgsTuple args;
		jam_context context { frame };
		jam_marshal_args<ArgsTuple, 0>::convert(args, context);
		// Invoke call, and return marshaled result.
		invoke<Return>(cxx_call,
			std::tuple_cat(
				std::make_tuple(jam_cxx_self::get<Class>(frame)), args),
			make_index_sequence<1 + std::tuple_size<ArgsTuple>::value> {});
		return L0;
	}
	catch (const std::exception &)
	{
		return L0;
	}
}

// Bound plain function...
template <class Call,
	class... Args,
	class Return,
	typename std::enable_if<!std::is_void<Return>::value, int>::type = 0>
static LIST * jam_call_function(
	FRAME * frame, int flags, Call cxx_call, Return (*)(Args...))
{
	using namespace mp;
	try
	{
		// Marshal arguments from frame->args.
		using ArgsTuple = std::tuple<typename remove_cvref<Args>::type...>;
		ArgsTuple args;
		jam_context context { frame };
		jam_marshal_args<ArgsTuple, 0>::convert(args, context);
		// Invoke call, and return marshaled result.
		return jam_binder::convert_to_bind_value<LIST *>(
			invoke<Return>(cxx_call, args,
				make_index_sequence<std::tuple_size<ArgsTuple>::value> {}));
	}
	catch (const std::exception &)
	{
		return L0;
	}
}
template <class Call,
	class... Args,
	class Return,
	typename std::enable_if<std::is_void<Return>::value, int>::type = 0>
static LIST * jam_call_function(
	FRAME * frame, int flags, Call cxx_call, Return (*)(Args...))
{
	using namespace mp;
	try
	{
		// Marshal arguments from frame->args.
		using ArgsTuple = std::tuple<typename remove_cvref<Args>::type...>;
		ArgsTuple args;
		jam_context context { frame };
		jam_marshal_args<ArgsTuple, 0>::convert(args, context);
		// Invoke call, and return marshaled result.
		invoke<Return>(cxx_call, args,
			make_index_sequence<std::tuple_size<ArgsTuple>::value> {});
		return L0;
	}
	catch (const std::exception &)
	{
		return L0;
	}
}

template <class... A>
struct jam_arg_spec_count_sum
{};
template <>
struct jam_arg_spec_count_sum<>
{
	// static constexpr std::size_t value = 0;
	enum : std::size_t
	{
		value = 0
	};
};
template <class X, class... A>
struct jam_arg_spec_count_sum<X, A...>
{
	// static constexpr std::size_t value
	//     = X::count + jam_arg_spec_count_sum<A...>::value;
	enum : std::size_t
	{
		value = X::count + jam_arg_spec_count_sum<A...>::value
	};
};

template <class... A>
struct jam_arg_spec_max_size
{
	// static constexpr std::size_t value
	//     = ::b2::bind::args_<A...>::count
	//     + jam_arg_spec_count_sum<A...>::value*2
	//     + 2;
	enum : std::size_t
	{
		value = ::b2::bind::args_<A...>::count
			+ jam_arg_spec_count_sum<A...>::value * 2 + 2
	};
};

// Builds the jam argument specifier for a rule from the bind arguments
// specification.
template <class... A>
struct jam_arg_spec
{
	using args_t = ::b2::bind::args_<A...>;

	// The jam specification, terminated by a `nullptr`. Note, the `spec` will
	// be larger than the actual specification. Hence use the `count` if one
	// needs to know the range.
	const char * spec[jam_arg_spec_max_size<A...>::value];

	// The number of strings in the specification, exclusive of the terminating
	// n`nullptr`.
	int count;

	jam_arg_spec(const args_t & args)
		: jam_arg_spec(args, mp::make_index_sequence<sizeof...(A)> {})
	{}

	template <std::size_t... I>
	jam_arg_spec(const args_t & args, mp::index_sequence<I...>)
	{
		count = 0;
		using std::get;
#if 0
        bool _[]{
            false,
            append_arg_list_to_spec(get<I>(args.arg))... };
#endif
		if (count == 0) spec[count++] = "*";
		spec[count] = nullptr;
	}

	template <int C>
	bool append_arg_list_to_spec(const ::b2::bind::arg_<C> & arg)
	{
		if (count > 0) spec[count++] = ":";
		for (int i = 0; i < C; ++i)
		{
			if (arg.args[i].count == ::b2::bind::param_::rest)
			{
				spec[count++] = "*";
				break;
			}
			spec[count++] = arg.args[i].name;
			switch (arg.args[i].count)
			{
				case ::b2::bind::param_::any: spec[count++] = "*"; break;
				case ::b2::bind::param_::many: spec[count++] = "+"; break;
				case ::b2::bind::param_::optional: spec[count++] = "?"; break;
			}
		}
		return true;
	}
};

template <class Return, class ArgSpec, class... Args>
void jam_native_bind(const string_t & module_name,
	const string_t & rule_name,
	ArgSpec & arg_spec,
	function_builtin_t native_rule,
	Return (*)(Args...))
{
	value_ptr rule_name_obj = object_new(rule_name.c_str());
	module_t * module = mod(module_name);

	int found = 0;
	RULE * const rule
		= (RULE *)hash_insert(demand_rules(module), rule_name_obj, &found);
	rule->name = object_copy(rule_name_obj);
	rule->procedure = 0;
	rule->module = module;
	rule->actions = 0;
	rule->exported = 0;
	// Register the function as a native jam rule.
	declare_native_rule(mod_cstr(module_name), rule_name.c_str(), arg_spec.spec,
		native_rule, 0);
	// Note, we don't check results of not finding the existing native rule,
	// for the obvious reason that we just created it a few lines above.
	native_rule_t * np
		= (native_rule_t *)hash_find(module->native_rules, rule_name_obj);
	// Set a global name for the native rule.
	function_set_rulename(
		np->procedure, object_new((module_name + "." + rule_name).c_str()));
	// Define the native rule in the class module.
	new_rule_body(module, np->name, np->procedure, 1);

	object_free(rule_name_obj);
}

template <class Return, class Class, class ArgSpec, class... Args>
void jam_bind(const string_t & module_name,
	const string_t & rule_name,
	ArgSpec & arg_spec,
	Return (Class::*method)(Args...))
{
	// out_printf( "jam_bind: %s.%s.\n", module_name.c_str(), rule_name.c_str()
	// );
	jam_native_bind(
		module_name, rule_name, arg_spec,
		[method](FRAME * frame, int flags) -> LIST * {
			return jam_call_method(
				frame, flags,
				[method](Class * self, Args... args) -> Return {
					return (self->*method)(args...);
				},
				static_cast<Return (*)(Class *, Args...)>(nullptr));
		},
		static_cast<Return (*)(Class *, Args...)>(nullptr));
}

template <class Return, class Class, class ArgSpec, class... Args>
void jam_bind(const string_t & module_name,
	const string_t & rule_name,
	ArgSpec & arg_spec,
	Return (Class::*method)(Args...) const)
{
	// out_printf( "jam_bind: %s.%s.\n", module_name.c_str(), rule_name.c_str()
	// );
	jam_native_bind(
		module_name, rule_name, arg_spec,
		[method](FRAME * frame, int flags) -> LIST * {
			return jam_call_method(
				frame, flags,
				[method](const Class * self, Args... args) -> Return {
					return (self->*method)(args...);
				},
				static_cast<Return (*)(Class *, Args...)>(nullptr));
		},
		static_cast<Return (*)(Class *, Args...)>(nullptr));
}

template <class Return, class ArgSpec, class... Args>
void jam_bind(const string_t & module_name,
	const string_t & rule_name,
	ArgSpec & arg_spec,
	Return (*function)(Args...))
{
	jam_native_bind(
		module_name, rule_name, arg_spec,
		[function](FRAME * frame, int flags) -> LIST * {
			return jam_call_function(
				frame, flags,
				[function](
					Args... args) -> Return { return (*function)(args...); },
				static_cast<Return (*)(Args...)>(nullptr));
		},
		static_cast<Return (*)(Args...)>(nullptr));
}

template <class Class, class ArgSpec, class... Args>
void jam_bind(const string_t & module_name,
	const string_t & rule_name,
	ArgSpec & arg_spec,
	Class *,
	::b2::bind::init_<Args...>)
{
	// out_printf( "jam_bind: %s.%s.\n", module_name.c_str(), rule_name.c_str()
	// );
	jam_native_bind(
		module_name, rule_name, arg_spec,
		[module_name, rule_name](FRAME * frame, int flags) -> LIST * {
			return jam_call_init(
				frame, flags,
				[](Args... args) -> Class * { return new Class(args...); },
				static_cast<Class * (*)(Args...)>(nullptr));
		},
		(void (*)(Args...)) nullptr);
}

void jam_binder::bind_module(const char * module_name) { mod(module_name); }

template <class Class, class... Bases>
void jam_binder::bind_class(const char * module_name,
	const char * class_name,
	::b2::bind::type_<Class>,
	::b2::bind::type_<Bases>...)
{
	value_ptr class_name_obj = object_new(class_name);
	LIST * class_name_list = list_new(class_name_obj);
	value_ptr bases_name_obj[]
		= { object_new(::b2::bind::class_<Bases, jam_binder>::name())...,
			  nullptr };
	LIST * bases_name_list = L0;
	for (value_ptr base_name_obj : bases_name_obj)
	{
		if (base_name_obj)
			bases_name_list = list_push_back(bases_name_list, base_name_obj);
	}
	/* value_ptr  class_module = */ make_class_module(
		class_name_list, bases_name_list);
	object_free(class_name_obj);
	for (value_ptr base_name_obj : bases_name_obj)
	{
		if (base_name_obj) object_free(base_name_obj);
	}
	list_free(bases_name_list);
}

template <class Function, class... A>
void jam_binder::bind_method(const char * module_name,
	const char * class_name,
	const char * method_name,
	::b2::bind::args_<A...> args,
	Function f)
{
	jam_arg_spec<A...> arg_spec { args };
	jam_bind(string_t("class@") + class_name, method_name, arg_spec, f);
}

template <class Class, class Init, class... A>
void jam_binder::bind_init(const char * module_name,
	const char * class_name,
	Class * c,
	Init i,
	::b2::bind::args_<A...> args)
{
	jam_arg_spec<A...> arg_spec { args };
	jam_bind(string_t("class@") + class_name, "__init__", arg_spec, c, i);
}

template <class Function, class... A>
void jam_binder::bind_function(const char * module_name,
	const char * function_name,
	::b2::bind::args_<A...> args,
	Function f)
{
	jam_arg_spec<A...> arg_spec { args };
	jam_bind(mod_str(module_name), function_name, arg_spec, f);
}

void jam_binder::eval_data(const char * module_name, const char * data)
{
	b2::jam::module_scope scope(
		context_ref.get<jam_context>().frame, module_name);
	parse_buffer(scope.name(), data, scope.frame());
}

void jam_binder::set_loaded(const char * module_name)
{
	b2::jam::variable modules_loaded("modules", ".loaded");
	modules_loaded += module_name;
}

template <>
void bind_jam(FRAME * f)
{
	jam_context context { f };
	jam_binder binder;
	binder.context_ref = context;
	binder
		// These have to be done in dependency order so that the init code
		// each executes is valid.
		.bind(jam::modules::modules_module())
		.bind(jam::builtin::root_module())
		.bind(jam::klass::class_module())
		.bind(jam::errors::errors_module())
		.bind(paths::paths_module())
		.bind(regex_module())
		.bind(set_module())
		.bind(string_module())
		.bind(sysinfo_module())
		.bind(version_module())
		.bind(db_module())
		.bind(command_db_module());
}

}} // namespace b2::jam
