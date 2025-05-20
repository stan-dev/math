// Copyright (c) 2022 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/cobalt/detail/util.hpp>

static_assert(boost::cobalt::detail::variadic_first<int, double, int>() == 1u);
static_assert(boost::cobalt::detail::variadic_first<const int &, int, double>() == std::numeric_limits<std::size_t>::max());
static_assert(boost::cobalt::detail::variadic_first<int, double>() == std::numeric_limits<std::size_t>::max());

static_assert(boost::cobalt::detail::get_variadic<0>(4.2, 3) == 4.2);
static_assert(boost::cobalt::detail::get_variadic<1>(4.2, 3) == 3);
static_assert(boost::cobalt::detail::get_variadic<0>(4, 2.) == 4u);
static_assert(boost::cobalt::detail::get_variadic<1>(4, 2.3) == 2.3);

static_assert(boost::cobalt::detail::variadic_has<int, double, int>);
static_assert(boost::cobalt::detail::variadic_has<int, int, double>);
static_assert(!boost::cobalt::detail::variadic_has<int, double>);

static_assert(std::is_same_v<boost::cobalt::detail::variadic_element_t<0u, double, int>, double>);
static_assert(std::is_same_v<boost::cobalt::detail::variadic_element_t<1u, int, int, double>, int>);
static_assert(std::is_same_v<boost::cobalt::detail::variadic_element_t<1u, int, double>, double>);
