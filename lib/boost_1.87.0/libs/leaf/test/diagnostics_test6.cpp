// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef BOOST_LEAF_TEST_SINGLE_HEADER
#   include "leaf.hpp"
#else
#   include <boost/leaf/config.hpp>
#   include <boost/leaf/diagnostics.hpp>
#   include <boost/leaf/on_error.hpp>
#   include <boost/leaf/result.hpp>
#endif

#if BOOST_LEAF_CFG_STD_STRING
#   include <sstream>
#   include <iostream>
#endif

#include <vector>

#include "lightweight_test.hpp"

namespace leaf = boost::leaf;

template <int> struct info { int value; };

leaf::result<void> f1()
{
    return leaf::new_error(info<1>{11}, info<2>{21});
}

leaf::result<void> f2()
{
    return leaf::try_handle_some(
        []
        {
            return f1();
        },
        []
        {
            return leaf::new_error(info<1>{12});
        } );
}

leaf::result<void> f3()
{
    return leaf::try_handle_some(
        []
        {
            return f1();
        },
        [](leaf::diagnostic_details const &)
        {
            return leaf::new_error(info<1>{13});
        } );
}

leaf::result<void> f4()
{
    return leaf::try_handle_some(
        []
        {
            return f1();
        },
        [](leaf::diagnostic_details const & e)
        {
            return e.error().load(info<1>{14});
        } );
}

int main()
{
    leaf::try_handle_all([]() -> leaf::result<void>
    {
        return f2();
    },
    [](leaf::diagnostic_details const & e)
    {
#if BOOST_LEAF_CFG_STD_STRING
        std::ostringstream st;
        st << e;
        std::string s = st.str();
        std::cout << s << std::endl;
        if( BOOST_LEAF_CFG_DIAGNOSTICS && BOOST_LEAF_CFG_CAPTURE )
        {
            BOOST_TEST_EQ(s.find(": 11"), s.npos);
            BOOST_TEST_NE(s.find(": 12"), s.npos);
            BOOST_TEST_EQ(s.find(": 21"), s.npos);
        }
#endif
    } );

    leaf::try_handle_all([]() -> leaf::result<void>
    {
        return f3();
    },
    [](leaf::diagnostic_details const & e)
    {
#if BOOST_LEAF_CFG_STD_STRING
        std::ostringstream st;
        st << e;
        std::string s = st.str();
        std::cout << s << std::endl;
        if( BOOST_LEAF_CFG_DIAGNOSTICS && BOOST_LEAF_CFG_CAPTURE )
        {
            BOOST_TEST_EQ(s.find(": 11"), s.npos);
            BOOST_TEST_NE(s.find(": 13"), s.npos);
            BOOST_TEST_EQ(s.find(": 21"), s.npos);
        }
#endif
    } );

    leaf::try_handle_all([]() -> leaf::result<void>
    {
        return f4();
    },
    [](leaf::diagnostic_details const & e)
    {
#if BOOST_LEAF_CFG_STD_STRING
        std::ostringstream st;
        st << e;
        std::string s = st.str();
        std::cout << s << std::endl;
        if( BOOST_LEAF_CFG_DIAGNOSTICS && BOOST_LEAF_CFG_CAPTURE )
        {
            BOOST_TEST_EQ(s.find(": 12"), s.npos);
            BOOST_TEST_NE(s.find(": 14"), s.npos);
            BOOST_TEST_NE(s.find(": 21"), s.npos);
        }
#endif
    } );

    return boost::report_errors();
}
