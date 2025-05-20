/*
Copyright 2017 Glen Joseph Fernandes
(glenjofe@gmail.com)

Distributed under the Boost Software License, Version 1.0.
(http://www.boost.org/LICENSE_1_0.txt)
*/
#include <boost/smart_ptr/make_local_shared.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>

int main()
{
    {
        boost::local_shared_ptr<int[]> result =
            boost::make_local_shared<int[]>(4, 1);
        BOOST_TEST(result[0] == 1);
        BOOST_TEST(result[1] == 1);
        BOOST_TEST(result[2] == 1);
        BOOST_TEST(result[3] == 1);
    }
    {
        boost::local_shared_ptr<int[4]> result =
            boost::make_local_shared<int[4]>(1);
        BOOST_TEST(result[0] == 1);
        BOOST_TEST(result[1] == 1);
        BOOST_TEST(result[2] == 1);
        BOOST_TEST(result[3] == 1);
    }
    {
        boost::local_shared_ptr<const int[]> result =
            boost::make_local_shared<const int[]>(4, 1);
        BOOST_TEST(result[0] == 1);
        BOOST_TEST(result[1] == 1);
        BOOST_TEST(result[2] == 1);
        BOOST_TEST(result[3] == 1);
    }
    {
        boost::local_shared_ptr<const int[4]> result =
            boost::make_local_shared<const int[4]>(1);
        BOOST_TEST(result[0] == 1);
        BOOST_TEST(result[1] == 1);
        BOOST_TEST(result[2] == 1);
        BOOST_TEST(result[3] == 1);
    }
    return boost::report_errors();
}
