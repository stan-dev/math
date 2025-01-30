/*
Copyright 2022-2023 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef BOOST_JAM_VALUE_H
#define BOOST_JAM_VALUE_H

#include "config.h"

#include "strview.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <string>

namespace b2 {

struct object
{
	virtual ~object() {}
};

struct value
{
	struct str_view
	{
		const char * str;
		std::size_t size;
		inline const char * begin() const { return str; }
		inline const char * end() const { return str + size; }
	};

	enum class type : char
	{
		null = '\0',
		string = 's',
		number = 'n',
		object = 'o'
	};

	mutable std::uint64_t hash64 = 0;

	virtual ~value() {}
	virtual type get_type() const = 0;
	virtual bool equal_to(const value & o) const = 0;
	virtual int compare_to(const value & o) const = 0;
	virtual str_view as_string() const = 0;
	virtual double as_number() const = 0;
	virtual object * as_object() const = 0;
	virtual value * to_string() = 0;
	inline const char * str() { return as_string().str; }
	static value * make(const char * str, std::size_t size);
	static inline value * make(const char * str)
	{
		return make(str, str ? strlen(str) : 0);
	}
	static inline value * make(const char * begin, const char * end)
	{
		return make(begin, end - begin);
	}
	static inline value * make(const std::string & str)
	{
		return make(str.c_str(), str.size());
	}
	static inline value * make(string_view str)
	{
		return make(str.data(), str.length());
	}
	static value * make(object * obj);
	static value * make(double v);
	static inline value * copy(value * v) { return v; }
	static inline void free(value *& v) { v = nullptr; }
	static void done();

	inline bool has_value() const { return get_type() == type::null; }

	template <typename T>
	static inline value * as_string(T v)
	{
		return value::make(std::to_string(v));
	}

	template <typename... Args>
	static inline value * format(const char * f, Args... args)
	{
		auto size = std::snprintf(nullptr, 0, f, args...);
		std::unique_ptr<char[]> s(new char[size + 1]);
		std::snprintf(s.get(), size + 1, f, args...);
		return value::make(s.get(), size);
	}

	inline int compare(const value & other) const
	{
		if (get_type() != other.get_type())
			return int(get_type()) - int(other.get_type());
		return compare_to(other);
	}
};

typedef value * value_ptr;

struct value_ref
{
	inline value_ref()
		: val(nullptr)
	{}
	inline value_ref(const value_ref & a)
		: val(value::copy(a.val))
	{}
	inline value_ref(value_ptr a)
		: val(value::copy(a))
	{}
	inline value_ref(const char * s)
		: val(value::make(s))
	{}
	inline value_ref(const std::string & s)
		: val(value::make(s.c_str()))
	{}

	inline ~value_ref()
	{
		if (val) value::free(val);
	}
	inline value_ptr release()
	{
		value_ptr r = val;
		val = nullptr;
		return r;
	}

	inline operator value_ptr() const { return val; }
	inline operator std::string() const { return val->str(); }
	inline bool operator==(value_ptr b) const { return val->equal_to(*b); }
	inline bool operator==(const char * v) const
	{
		return (*this) == value_ref(v);
	}
	template <typename T>
	inline bool operator!=(T v) const
	{
		return !((*this) == v);
	}
	inline value_ptr operator->() const { return val; }
	inline value & operator*() const { return *val; }
	inline bool has_value() const { return val != nullptr; }

	inline value_ref operator+(const std::string & v) const
	{
		return value_ref((*this)->str() + v);
	}

	friend inline value_ref operator+(const std::string & a, value_ref b)
	{
		const char * c = b.has_value() ? b->str() : "";
		return value_ref(a + c);
	}

	struct hash_function
	{
		inline std::size_t operator()(const value_ref & a) const
		{
			return (std::size_t)(
				(std::numeric_limits<std::size_t>::max)() & a->hash64);
		}
	};

	struct equal_function
	{
		inline bool operator()(const value_ref & a, const value_ref & b) const
		{
			return (a->hash64 == b->hash64) && a->equal_to(*b);
		}
	};

	private:
	value_ptr val = nullptr;
};

} // namespace b2

#endif
