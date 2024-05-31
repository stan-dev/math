//
// Copyright (c) 2023 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <concepts>
static_assert(!std::derived_from<int, double>);
static_assert(std::same_as<int, int>);
static_assert(std::convertible_to<int, double>);
