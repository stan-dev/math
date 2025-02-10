// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/leaf/config.hpp>

#if !BOOST_LEAF_CFG_CAPTURE

#include <iostream>

int main()
{
    std::cout << "Unit test not applicable." << std::endl;
    return 0;
}

#else

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

int count = 0;

template <int N>
struct info
{
    info() noexcept
    {
        ++count;
    }

    info( info const & ) noexcept
    {
        ++count;
    }

    ~info() noexcept
    {
        --count;
    }

    template <class CharT, class Traits>
    friend std::ostream & operator<<( std::basic_ostream<CharT, Traits> & os, info const & )
    {
        return os << "info<" << N << "> instance";
    }
};

int main()
{
    {
        leaf::result<int> r = leaf::try_capture_all(
            []() -> leaf::result<int>
            {
                return leaf::new_error( info<1>{}, info<3>{} );
            } );
        BOOST_TEST_EQ(count, 2);

#if BOOST_LEAF_CFG_STD_STRING
        {
            std::ostringstream st;
            st << r;
            std::string s = st.str();
            std::cout << s << std::endl;
            if( BOOST_LEAF_CFG_DIAGNOSTICS )
            {
                BOOST_TEST_NE(s.find("info<1> instance"), s.npos);
                BOOST_TEST_NE(s.find("info<3> instance"), s.npos);
            }
        }
#endif

        int ret = leaf::try_handle_all(
            [&]
            {
                return std::move(r);
            },
            []( info<1>, info<3> )
            {
                return 42;
            },
            []
            {
                return -42;
            } );
        BOOST_TEST_EQ(ret, 42);
    }
    BOOST_TEST_EQ(count, 0);

    return boost::report_errors();
}

#endif
