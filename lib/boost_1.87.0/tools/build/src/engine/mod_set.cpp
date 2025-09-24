/*
Copyright 2019-2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "mod_set.h"

#include <algorithm>
#include <unordered_map>

void b2::set::add(list_cref values)
{
	for (auto v : values) elements.insert(value_ref(v));
}

void b2::set::add(const set & value)
{
	elements.insert(value.elements.cbegin(), value.elements.cend());
}

bool b2::set::contains(value_ref value) const
{
	return elements.find(value) != elements.end();
}

b2::list_ref b2::set::to_list() const
{
	list_ref result;
	for (auto v : elements) result.push_back(v);
	return result;
}

b2::list_ref b2::set::difference(b2::list_cref a, b2::list_cref b)
{
	// Somewhat less efficient variation on difference to maintain the order
	// of sequence "a".
	value_set rhs { b.begin(), b.end() };
	list_ref result;
	for (auto val : a)
		if (rhs.count(val) == 0) result.push_back(val);
	return result;
}

b2::list_ref b2::set::intersection(b2::list_cref a, b2::list_cref b)
{
	value_set lhs { a.begin(), a.end() };
	list_ref result;
	for (auto val : b)
		if (lhs.count(val) > 0) result.push_back(val);
	return result;
}

bool b2::set::equal(b2::list_cref a, b2::list_cref b)
{
	std::unordered_map<value_ref, unsigned, value_ref::hash_function,
		value_ref::equal_function>
		val_count;
	for (auto val : a) val_count[val] += 1;
	for (auto val : b) val_count[val] += 1;
	for (auto & vc : val_count)
		if (vc.second == 1) return false;
	return true;
}

const char * b2::set_module::init_code = R"jam(
rule __test__ ( )
{
    import assert ;

    assert.result 0 1 4 6 8 9 : difference   0 1 2 3 4 5 6 7 8 9 : 2 3 5 7 ;
    assert.result 2 5 7       : intersection 0 1 2   4 5 6 7 8 9 : 2 3 5 7 ;

    assert.true  equal         :         ;
    assert.true  equal 1 1 2 3 : 3 2 2 1 ;
    assert.false equal 2 3     : 3 2 2 1 ;
}
)jam";
