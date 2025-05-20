/*=============================================================================
    Copyright (c) 2024 Nikita Kniazev

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#include <boost/spirit/home/x3.hpp>

#define OP(x) \
template <typename P> \
auto adl_checker_impl(P& p, int) -> decltype(void(x)) { \
    static_assert(sizeof(P) == 0, #x " is unconstrained"); \
}

template <typename T>
auto adl_checker_impl(T&, long) {}

OP(p >> p)
#if defined(_MSC_VER) && _MSC_VER >= 1910
OP(p > p)
OP(p | p)
OP(p - p)
OP(p % p)
OP(-p)
OP(!p)
OP(*p)
OP(+p)
//OP(&p)
#endif

template <typename P>
void adl_checker(P& p)
{
    adl_checker_impl(p, 0);
}

namespace boost { namespace spirit { namespace x3 {
    struct test_base {};
}}}

struct not_a_parser : boost::spirit::x3::test_base {};

int main()
{
    not_a_parser test;
    adl_checker(test);
}
