/*
Copyright 2019-2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "mod_version.h"
#include "modules.h"
#include <algorithm>

namespace b2 {

bool version_less(const std::vector<int> & lhs, const std::vector<int> & rhs)
{
	auto left = lhs.begin();
	auto left_end = lhs.end();
	auto right = rhs.begin();
	auto right_end = rhs.end();
	while (left != left_end || right != right_end)
	{
		int left_value = 0;
		int right_value = 0;
		if (left != left_end) left_value = *(left++);
		if (right != right_end) right_value = *(right++);
		if (left_value > right_value) return false;
		if (left_value < right_value) return true;
	}
	return false;
}

} // namespace b2
