/**
 *   Copyright (C) 2018 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/parser/parser.hpp>

#include <boost/core/lightweight_test.hpp>


using namespace boost::parser;


auto true_ = [](auto & context) { return true; };
auto false_ = [](auto & context) { return false; };

auto three = [](auto & context) { return 3; };


int main()
{

// if_
{
    {
        std::string str = "a";
        auto result = parse(str, if_(true_)[char_]);
        BOOST_TEST(result);
        BOOST_TEST(*result == 'a');
        BOOST_TEST(!parse(str, if_(false_)[char_]));
    }
    {
        std::string str = "a";
        auto result = parse(str, if_(true_)[char_]);
        BOOST_TEST(result);
        BOOST_TEST(*result == 'a');
        BOOST_TEST(!parse(str, if_(false_)[char_]));
        BOOST_TEST(!parse(str, if_(true_)[char_('b')]));
    }
}

// switch_
{
    {
        std::string str = "a";
        auto result = parse(str, switch_(true_)(true, char_));
        BOOST_TEST(result);
        BOOST_TEST(*result == 'a');
    }
    {
        std::string str = "a";
        auto result = parse(str, switch_(true)(true_, char_));
        BOOST_TEST(result);
        BOOST_TEST(*result == 'a');
    }
    {
        std::string str = "a";
        auto result = parse(str, switch_(true)(true, char_));
        BOOST_TEST(result);
        BOOST_TEST(*result == 'a');
    }
    {
        std::string str = "a";
        auto result = parse(str, switch_(true_)(true_, char_));
        BOOST_TEST(result);
        BOOST_TEST(*result == 'a');
    }

    {
        std::string str = "a";
        auto result = parse(str, switch_(true_)(false, char_));
        BOOST_TEST(!result);
    }
    {
        std::string str = "a";
        auto result = parse(str, switch_(true)(false_, char_));
        BOOST_TEST(!result);
    }
    {
        std::string str = "a";
        auto result = parse(str, switch_(true)(false, char_));
        BOOST_TEST(!result);
    }
    {
        std::string str = "a";
        auto result = parse(str, switch_(true_)(false_, char_));
        BOOST_TEST(!result);
    }

    {
        std::string str = "a";
        auto result =
            parse(str, switch_(three)(1, char_('b'))(three, char_('a')));
        BOOST_TEST(result);
        BOOST_TEST(*result == 'a');
    }
    {
        std::string str = "a";
        auto result = parse(str, switch_(3)(1, char_('b'))(three, char_('a')));
        BOOST_TEST(result);
        BOOST_TEST(*result == 'a');
    }
    {
        std::string str = "a";
        auto result = parse(str, switch_(3)(1, char_('b'))(3, char_('a')));
        BOOST_TEST(result);
        BOOST_TEST(*result == 'a');
    }
    {
        std::string str = "a";
        auto result = parse(str, switch_(three)(1, char_('b'))(3, char_('a')));
        BOOST_TEST(result);
        BOOST_TEST(*result == 'a');
    }

    {
        std::string str = "a";
        auto result =
            parse(str, switch_(three)(1, char_('a'))(three, char_('b')));
        BOOST_TEST(!result);
    }
    {
        std::string str = "a";
        auto result = parse(str, switch_(3)(1, char_('a'))(three, char_('b')));
        BOOST_TEST(!result);
    }
    {
        std::string str = "a";
        auto result = parse(str, switch_(3)(1, char_('a'))(3, char_('b')));
        BOOST_TEST(!result);
    }
    {
        std::string str = "a";
        auto result = parse(str, switch_(three)(1, char_('a'))(3, char_('b')));
        BOOST_TEST(!result);
    }
}

return boost::report_errors();
}
