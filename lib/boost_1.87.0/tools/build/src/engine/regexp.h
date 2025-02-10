/*
 * Definitions etc. for regexp(3) routines.
 *
 * Caveat:  this is V8 regexp(3) [actually, a reimplementation thereof],
 * not the System V one.
 */
#ifndef B2_REGEXP_H
#define B2_REGEXP_H

#include "config.h"

#include "mem.h"
#include "strview.h"
#include "types.h"
#include "value.h"

#include <functional>
#include <memory>

#define NSUBEXP 10

namespace b2 { namespace regex {

// The resulting matches for a regex match. Expression 0 is the full match.
// And expressions [1,NSUBEXP] are the subexpressions matched.
struct regex_expr
{
	string_view sub[NSUBEXP];
};

// The compiled regex program to match with.
struct regex_prog;

struct program
{
	program() = default;
	program(program && o)
		: compiled(std::move(o.compiled))
	{}
	program(const program &) = default;
	explicit program(const char * pattern);

	// types
	struct result_iterator;

	void reset(const char * pattern);

	result_iterator search(const string_view & str);
	result_iterator search(const char * str_begin, const char * str_end);
	result_iterator search(const char * str_begin);

	private:
	const regex_prog * compiled = nullptr;

	static regex_prog & compile(const char * patter);
};

struct program::result_iterator
{
	using iterator_category = std::forward_iterator_tag;
	using value_type = string_view;
	using difference_type = std::ptrdiff_t;
	using pointer = const value_type *;
	using reference = const value_type &;

	result_iterator(const regex_prog & c, const string_view & s);
	result_iterator(const result_iterator & o) = default;
	result_iterator(result_iterator && o)
		: compiled(std::move(o.compiled))
		, expressions(std::move(o.expressions))
		, rest(std::move(o.rest))
	{}

	inline result_iterator & operator++()
	{
		advance();
		return *this;
	}
	inline reference operator*() const { return (*this)[0]; }
	inline pointer operator->() const { return &(*this)[0]; }
	explicit inline operator bool() const { return !(*this)[0].empty(); }
	inline reference operator[](std::size_t i) const
	{
		static const value_type invalid { nullptr, 0 };
		return i <= NSUBEXP ? expressions.sub[i] : invalid;
	}

	private:
	const regex_prog * compiled = nullptr;
	regex_expr expressions;
	value_type rest;

	void advance();
};

inline program::result_iterator program::search(const string_view & str)
{
	return result_iterator(*compiled, str);
}

inline program::result_iterator program::search(
	const char * str_begin, const char * str_end)
{
	return this->search(string_view(str_begin, str_end));
}

inline program::result_iterator program::search(const char * str_begin)
{
	return this->search(string_view(str_begin, std::strlen(str_begin)));
}

}} // namespace b2::regex

#endif
