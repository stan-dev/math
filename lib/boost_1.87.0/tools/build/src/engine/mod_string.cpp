/*
Copyright 2019-2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "mod_string.h"

#include <algorithm>

static const char * whitespace_chars = " \t\n";

b2::value_ref b2::string_whitespace()
{
	return b2::value_ref(b2::value::make(whitespace_chars));
}

b2::list_ref b2::string_chars(b2::value_ref string)
{
	list_ref chars;
	for (auto ch : string->as_string()) chars.push_back(value::make(&ch, 1));
	return chars;
}

b2::value_ref b2::string_abbreviate(b2::value_ref string)
{
	// Anything less than 4 characters gets no abbreviation.
	if (string->as_string().size < 4) return string;
	std::string s = string->str();
	// Separate the initial letter in case it's a vowel.
	auto first_c = s[0];
	s.erase(0, 1);
	// Drop trailing "ing".
	if (s.compare(s.size() - 3, 3, "ing") == 0) s.erase(s.size() - 3);
	// Reduce all doubled characters to one.
	s.erase(std::unique(s.begin(), s.end()), s.end());
	// Chop all vowels out of the remainder.
	s.erase(std::remove_if(s.begin(), s.end(),
				[](std::string::value_type c) {
					// Sorted by English vowel frequency.
					return c == 'e' || c == 'E' || c == 'a' || c == 'A'
						|| c == 'o' || c == 'O' || c == 'i' || c == 'I'
						|| c == 'u' || c == 'U';
				}),
		s.end());
	// Shorten remaining consonants to 4 characters.
	if (s.size() > 4) s.erase(4);
	// Glue the initial character back on to the front.
	s.insert(0, 1, first_c);

	return value_ref(value::make(s.c_str()));
}

b2::value_ref b2::string_join(b2::list_cref strings, b2::value_ref separator)
{
	if (strings.empty()) return value_ref(value::make(""));
	std::string result = (*strings.begin())->str();
	auto strings_i = ++strings.begin();
	auto strings_e = strings.end();
	for (; strings_i != strings_e; ++strings_i)
	{
		if (separator.has_value()) result += separator->str();
		result += (*strings_i)->str();
	}
	return value_ref(value::make(result.c_str()));
}

b2::list_ref b2::string_words(std::string s, b2::list_cref whitespace)
{
	std::string delimiters { whitespace_chars };
	if (!whitespace.empty())
		delimiters = string_join(whitespace, value_ref(value::make("")));
	list_ref words;
	std::string::size_type word_begin = 0;
	std::string::size_type word_ws = s.find_first_of(delimiters);
	while (word_ws != std::string::npos)
	{
		words.push_back(s.substr(word_begin, word_ws - word_begin));
		word_begin = s.find_first_not_of(delimiters, word_ws + 1);
		word_ws = s.find_first_of(delimiters, word_begin);
	}
	if (word_begin < s.size()) words.push_back(s.substr(word_begin));
	return words;
}

bool b2::string_is_whitespace(value_ref s)
{
	if (!s.has_value()) return true;
	if (s->get_type() != value::type::string) return false;
	for (auto c : s->as_string())
	{
		if (c != ' ' && c != '\t' && c != '\n') return false;
	}
	return true;
}

const char * b2::string_module::init_code = R"jam(
rule __test__ ( )
{
    import assert ;
    assert.result a b c : chars abc ;

    assert.result rntm : abbreviate runtime ;
    assert.result ovrld : abbreviate overload ;
    assert.result dbg : abbreviate debugging ;
    assert.result async : abbreviate asynchronous ;
    assert.result pop : abbreviate pop ;
    assert.result aaa : abbreviate aaa ;
    assert.result qck : abbreviate quack ;
    assert.result sttc : abbreviate static ;

    # Check boundary cases.
    assert.result a : chars a ;
    assert.result : chars "" ;
    assert.result a b c d e f g h : chars abcdefgh ;
    assert.result a b c d e f g h i : chars abcdefghi ;
    assert.result a b c d e f g h i j : chars abcdefghij ;
    assert.result a b c d e f g h i j k : chars abcdefghijk ;

    assert.result a//b/c/d : join a "" b c d : / ;
    assert.result abcd : join  a "" b c d ;

    assert.result a b c : words "a b	c" ;

    assert.true is-whitespace "     	" ;
    assert.false is-whitespace "  a b c	" ;
    assert.true is-whitespace "" ;
    assert.true is-whitespace ;
}
)jam";
