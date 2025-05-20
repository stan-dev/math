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

#include "lightweight_test.hpp"

namespace leaf = boost::leaf;

struct info { int value; };

leaf::result<void> f1()
{
    return leaf::new_error(info{41});
}

leaf::result<void> f2()
{
    auto load = leaf::on_error(info{42});
    return f1();
}

int main()
{
    leaf::try_handle_all(
        []() -> leaf::result<void>
        {
            return f2();
        },
        []( leaf::diagnostic_details const & di )
        {
#if BOOST_LEAF_CFG_STD_STRING
            std::ostringstream st;
            st << di;
            std::string s = st.str();
            std::cout << s << std::endl;
            if( BOOST_LEAF_CFG_DIAGNOSTICS && BOOST_LEAF_CFG_CAPTURE )
            {
                BOOST_TEST_NE(s.find("41"), s.npos);
                BOOST_TEST_EQ(s.find("42"), s.npos);
            }
#endif
        } );

    return boost::report_errors();
}
