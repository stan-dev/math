/*=============================================================================
    Copyright (c) 2005-2007 Joel de Guzman
    Copyright (c) 2021 Beojan Stanislaus

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#include <iostream>
#if __cplusplus >= 201402L || (defined(_MSC_FULL_VER) && _MSC_FULL_VER >= 190024210)

#include <boost/detail/lightweight_test.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/stl/tuple.hpp>

#include <tuple>

namespace {
void tuple_get_test() {
    using boost::phoenix::get_;
    using boost::phoenix::arg_names::arg1;
    std::tuple<int, double, char> tup{1, 1.6, 'b'};
    BOOST_TEST(get_<0>(arg1)(tup) == std::get<0>(tup));
    BOOST_TEST(get_<1>(arg1)(tup) == std::get<1>(tup));
    BOOST_TEST(get_<2>(arg1)(tup) == std::get<2>(tup));
    BOOST_TEST(get_<int>(arg1)(tup) == std::get<int>(tup));
    BOOST_TEST(get_<double>(arg1)(tup) == std::get<double>(tup));
    BOOST_TEST(get_<char>(arg1)(tup) == std::get<char>(tup));
    return;
}
} // namespace

int main() {
    tuple_get_test();
    return boost::report_errors();
}

#else
int main() {
    std::cout << "Skipping test -- Tuple support requires C++14 or above" << std::endl;
    return 0;
}
#endif
