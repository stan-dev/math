//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/test/unit_test.hpp>

#ifdef BOOST_TEST_ALTERNATIVE_INIT_API
int main(int argc, char* argv[])
{
    return ::boost::unit_test::unit_test_main([] { return true; }, argc, argv);
}
#else
::boost::unit_test::test_suite* init_unit_test_suite(int, char*[]) { return nullptr; }
#endif
