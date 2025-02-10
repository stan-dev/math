// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef BOOST_LEAF_TEST_SINGLE_HEADER
#   include "leaf.hpp"
#else
#   include <boost/leaf/result.hpp>
#   include <boost/leaf/handle_errors.hpp>
#endif

#if BOOST_LEAF_CFG_STD_STRING
#   include <sstream>
#   include <iostream>
#endif

#include "lightweight_test.hpp"

namespace leaf = boost::leaf;

struct non_printable_value
{
};

struct e_err
{
    template <class CharT, class Traits>
    friend std::ostream & operator<<( std::basic_ostream<CharT, Traits> & os, e_err const & )
    {
        return os << "e_err";
    }
};

int main()
{
    {
        leaf::result<int> r = 42;
        BOOST_TEST(r);
#if BOOST_LEAF_CFG_STD_STRING
        std::ostringstream ss;
        ss << r;
        std::string s = ss.str();
        std::cout << s << std::endl;
        if( BOOST_LEAF_CFG_DIAGNOSTICS )
            BOOST_TEST_EQ(s, "42");
#endif
    }

    {
        leaf::result<non_printable_value> r;
        BOOST_TEST(r);
#if BOOST_LEAF_CFG_STD_STRING
        std::ostringstream ss;
        ss << r;
        std::string s = ss.str();
        std::cout << s << std::endl;
        if( BOOST_LEAF_CFG_DIAGNOSTICS )
            BOOST_TEST_EQ(s, "{not printable}");
#endif
    }

    {
        leaf::result<void> r;
        BOOST_TEST(r);
#if BOOST_LEAF_CFG_STD_STRING
        std::ostringstream ss;
        ss << r;
        std::string s = ss.str();
        std::cout << s << std::endl;
        if( BOOST_LEAF_CFG_DIAGNOSTICS )
            BOOST_TEST_EQ(s, "No error");
#endif
    }

    {
        leaf::result<int> r = leaf::new_error(e_err{ });
        BOOST_TEST(!r);
#if BOOST_LEAF_CFG_STD_STRING
        std::ostringstream ss;
        ss << r;
        std::string s = ss.str();
        std::cout << s << std::endl;
        leaf::error_id err = r.error();
        if( BOOST_LEAF_CFG_DIAGNOSTICS )
            BOOST_TEST_EQ(s, "Error serial #" + std::to_string(err.value()/4));
#endif
    }

#if BOOST_LEAF_CFG_CAPTURE
    {
        leaf::result<int> r = leaf::try_capture_all(
            []() -> leaf::result<int>
            {
                return leaf::new_error(e_err{ });
            } );
#if BOOST_LEAF_CFG_STD_STRING
        std::ostringstream ss;
        ss << r;
        std::string s = ss.str();
        std::cout << s << std::endl;
        leaf::error_id err = r.error();
        BOOST_TEST_NE(s.find("Error serial #" + std::to_string(err.value()/4)), s.npos);
        if( BOOST_LEAF_CFG_DIAGNOSTICS )
        {
            BOOST_TEST_NE(s.find("Captured:"), s.npos);
            BOOST_TEST_NE(s.find("e_err: e_err"), s.npos);
        }
#endif
    }
#endif

    return boost::report_errors();
}
