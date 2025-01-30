/*
Copyright 2022-2023 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "value.h"

#include "jam.h"
#include "mem.h"
#include "mp.h"
#include "object.h"
#include "output.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <unordered_set>
#include <vector>

static inline void hash_fnv1a(
	std::uint64_t & value, void * data, std::size_t size)
{
	value = 0xcbf29ce484222325;
	std::uint8_t * d = static_cast<std::uint8_t *>(data);
	std::uint8_t * e = d + size;
	while (d != e) value = (value * 0x100000001B3) ^ *(d++);
}
inline void hash_fnv1a(std::uint32_t & value, void * data, std::size_t size)
{
	value = 0x811c9dc5;
	std::uint8_t * d = static_cast<std::uint8_t *>(data);
	std::uint8_t * e = d + size;
	while (d != e) value = (value * 0x1000193) ^ *(d++);
}
#if 0
static inline void hash_fnv1a(std::uint32_t & value32,
	std::uint64_t & value64,
	void * data,
	std::size_t size)
{
	value32 = 0x811c9dc5;
	value64 = 0xcbf29ce484222325;
	std::uint8_t * d = static_cast<std::uint8_t *>(data);
	std::uint8_t * e = d + size;
	while (d != e)
	{
		value32 = (value32 * 0x1000193) ^ *d;
		value64 = (value64 * 0x100000001B3) ^ *d;
		d += 1;
	}
}
#endif
static inline bool str_view_equal(const char * str_a,
	std::size_t size_a,
	const char * str_b,
	std::size_t size_b)
{
	// Some MSVC versions with the debug iterators will complain about using
	// std::equal with raw pointers. So we fall back to what std::equal
	// itself calls, i.e. we call memcmp instead of std::equal. And this is
	// how good C++ dies.
	return (str_a == str_b)
		|| ((str_a != nullptr) && (str_b != nullptr) && (size_a == size_b)
			&& (std::memcmp(str_a, str_b, size_a) == 0));
}
static inline int str_view_cmp(const char * str_a,
	std::size_t size_a,
	const char * str_b,
	std::size_t size_b)
{
	if (str_a == str_b) return 0;
	if (str_a == nullptr) return -1;
	if (str_b == nullptr) return 1;
	int result = std::memcmp(str_a, str_b, size_a < size_b ? size_a : size_b);
	if (result == 0) return size_a < size_b ? -1 : 1;
	return result;
}

namespace b2 {

struct value_base : value
{
	virtual type get_type() const override { return type::null; }
	virtual str_view as_string() const override { return { nullptr, 0 }; }
	virtual double as_number() const override
	{
		using nl = std::numeric_limits<double>;
		return nl::has_quiet_NaN ? nl::quiet_NaN() : nl::lowest();
	}
	virtual object * as_object() const override { return nullptr; }
};

struct value_str : value_base
{
	using value_type = char[8];

	static value_str * make(const char * v, std::size_t n)
	{
		const auto value_size = sizeof(value_str::value_type);
		std::size_t m = n + 1 < value_size ? value_size : n + 1;
		void * o = BJAM_MALLOC(sizeof(value_str) - value_size + m);
		return b2::jam::ctor_ptr<value_str>(o, v, n);
	}
	value_str(const char * v, std::size_t n)
		: size(n)
	{
		char * vv = value;
		std::memcpy(vv, v, n);
		vv[n] = '\0';
		hash_fnv1a(hash64, vv, n);
	}
	virtual type get_type() const override { return type::string; }
	virtual bool equal_to(const value & o) const override
	{
		return str_view_equal(
			value, size, o.as_string().str, o.as_string().size);
	}
	virtual int compare_to(const value & o) const override
	{
		return str_view_cmp(value, size, o.as_string().str, o.as_string().size);
	}
	virtual str_view as_string() const override { return { value, size }; }
	virtual value * to_string() override { return this; }

	std::size_t size;
	value_type value;
};

struct value_number : value_base
{
	static value_number * make(double v)
	{
		return b2::jam::make_ptr<value_number>(v);
	}
	value_number(double v)
		: value(v)
	{
		hash_fnv1a(hash64, (void *)&value, sizeof(value));
	}
	virtual type get_type() const override { return type::number; }
	virtual bool equal_to(const value & o) const override
	{
		return as_number() == o.as_number();
	}
	virtual int compare_to(const value & o) const override
	{
		if (as_number() == o.as_number()) return 0;
		if (as_number() < o.as_number()) return -1;
		return 1;
	}
	virtual double as_number() const override { return value; }
	virtual value * to_string() override
	{
		return value::make(std::to_string(value));
	}
	double value;
};

struct value_object : value_base
{
	static value_object * make(object * v)
	{
		return b2::jam::make_ptr<value_object>(v);
	}
	value_object(object * v)
		: value(v)
	{
		hash_fnv1a(hash64, (void *)&v, sizeof(v));
	}
	virtual type get_type() const override { return type::object; }
	virtual bool equal_to(const value & o) const override
	{
		return as_object() == o.as_object();
	}
	virtual int compare_to(const value & o) const override
	{
		if (as_object() == o.as_object()) return 0;
		if (as_object() < o.as_object()) return -1;
		return 1;
	}
	virtual object * as_object() const override { return value.get(); }
	virtual value * to_string() override
	{
		return value::make("@" + std::to_string(std::uintptr_t(value.get())));
	}
	std::unique_ptr<object> value;
};

struct value_str_view : value_base
{
	value_str_view(const char * v, std::size_t n)
		: value(v)
		, size(n)
	{
		hash_fnv1a(hash64, (void *)v, n);
	}
	virtual ~value_str_view() {}
	virtual bool equal_to(const value & o) const override
	{
		return str_view_equal(
			value, size, o.as_string().str, o.as_string().size);
	}
	virtual int compare_to(const value & o) const override
	{
		return str_view_cmp(
			value, size, o.as_string().str, o.as_string().size);
	}
	virtual str_view as_string() const override { return { value, size }; }
	virtual value * to_string() override { return nullptr; }

	const char * value = nullptr;
	std::size_t size = 0;
};

struct value_hash_f
{
	inline std::size_t operator()(const value * o) const
	{
		return (std::size_t)(
			(std::numeric_limits<std::size_t>::max)() & o->hash64);
	}
};

struct value_eq_f
{
	inline bool operator()(const value * a, const value * b) const
	{
		return (a->hash64 == b->hash64) && a->equal_to(*b);
	}
};

struct safe_value_cache
{
	template <typename Test, typename Val, typename A0, typename... An>
	typename std::enable_if<
		!std::is_same<object *, typename mp::remove_cvref<A0>::type>::value,
		value_ptr>::type
		save(A0 a0, An... an)
	{
		Test test_val(a0, an...);
		std::lock_guard<std::mutex> guard(mutex);
		auto existing = cache.find(&test_val);
		if (existing != cache.end()) return *existing;
		value_ptr result = Val::make(a0, an...);
		cache.insert(result);
		return result;
	}

	template <typename Test, typename Val>
	value_ptr save(object * obj)
	{
		value_object test_val(obj);
		std::lock_guard<std::mutex> guard(mutex);
		auto existing = cache.find(&test_val);
		if (existing != cache.end()) return *existing;
		value_ptr result = value_object::make(test_val.value.release());
		cache.insert(result);
		return result;
	}

	void reset()
	{
		std::lock_guard<std::mutex> guard(mutex);
		for (value * o : cache)
		{
			b2::jam::free_ptr(o);
		}
		cache.clear();
	}

	private:
	using value_cache_t = std::unordered_set<value *, value_hash_f, value_eq_f>;
	value_cache_t cache;
	std::mutex mutex;
};

static safe_value_cache & value_cache()
{
	static safe_value_cache c;
	return c;
}

value_ptr value::make(const char * string, std::size_t size)
{
	if (!string || size == 0)
	{
		string = "";
		size = 0;
	}
	return value_cache().save<value_str_view, value_str>(string, size);
}

value_ptr value::make(object * obj)
{
	return value_cache().save<value_object, value_object>(obj);
}

value_ptr value::make(double v)
{
	return value_cache().save<value_number, value_number>(v);
}

void value::done(void) { value_cache().reset(); }

} // namespace b2
