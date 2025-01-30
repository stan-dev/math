// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef BOOST_LEAF_TEST_SINGLE_HEADER
#   include "leaf.hpp"
#else
#   include <boost/leaf/on_error.hpp>
#   include <boost/leaf/diagnostics.hpp>
#   include <boost/leaf/result.hpp>
#endif

#if BOOST_LEAF_CFG_STD_STRING
#   include <sstream>
#   include <iostream>
#endif

#include "lightweight_test.hpp"

namespace leaf = boost::leaf;

int cleanup_count;

int global;

int get_global() noexcept
{
    return global;
}

template <int>
struct info
{
    int value;
};

leaf::error_id g()
{
    global = 0;
    auto load = leaf::on_error( []{++cleanup_count;}, []{ return info<42>{get_global()}; }, []{ return info<-42>{-42}; } );
    global = 42;
    return leaf::new_error();
}

leaf::error_id f()
{
    auto load = leaf::on_error([]{++cleanup_count;});
    return g();
}

int main()
{
    int r = leaf::try_handle_all(
        []() -> leaf::result<int>
        {
            return f();
        },
        []( info<42> const & i42, leaf::diagnostic_info const & di )
        {
            BOOST_TEST_EQ(i42.value, 42);
#if BOOST_LEAF_CFG_STD_STRING
            std::ostringstream ss;
            ss << di;
            std::string s = ss.str();
            std::cout << s << std::endl;
            if( BOOST_LEAF_CFG_DIAGNOSTICS )
                BOOST_TEST_NE(s.find("info<42>"), s.npos);
#endif
            return 1;
        },
        []
        {
            return 2;
        } );
    BOOST_TEST_EQ(r, 1);
    BOOST_TEST_EQ(cleanup_count, 2);
    return boost::report_errors();
}
