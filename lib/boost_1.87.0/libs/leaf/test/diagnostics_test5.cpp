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

namespace
{
    int counter = 0;
}

template <int N>
struct info
{
    info(info const &) = delete;
    info & operator=(info const &) = delete;
    int acc = 0;
    info()
    {
        ++counter;
    }
    info( info && x ):
        acc(x.acc)
    {
        ++counter;
    }
    ~info()
    {
        --counter;
    }
    void accumulate()
    {
        ++acc;
    }
    friend std::ostream & operator<<( std::ostream & os, info const & x )
    {
        return os << "info<" << N << ">: acc=" << x.acc;
    }
};

leaf::result<void> f1()
{
    return leaf::new_error(info<1>(), [](info<4> & x){ x.accumulate(); });
}

leaf::result<void> f2()
{
    auto load = leaf::on_error(info<2>{}, [](){ return info<3>(); }, [](info<4> & x){ x.accumulate(); });
    return f1();
}

leaf::result<void> f3()
{
    return leaf::try_handle_some(
        []() -> leaf::result<void>
        {
            return f2();
        },
        []( leaf::diagnostic_details const & e )
        {
            return e.error();
        } );
}

int main()
{
    BOOST_TEST_EQ(counter, 0);
    leaf::try_handle_all(
        []() -> leaf::result<void>
        {
            return f3();
        },
        []( info<1> const &, leaf::diagnostic_details const & di )
        {
#if BOOST_LEAF_CFG_STD_STRING
            std::ostringstream st;
            st << di;
            std::string s = st.str();
            std::cout << s << std::endl;
            if( BOOST_LEAF_CFG_DIAGNOSTICS && BOOST_LEAF_CFG_CAPTURE )
            {
                auto const n1 = s.find("info<1>: acc=0");
                auto const n2 = s.find("info<2>: acc=0");
                auto const n3 = s.find("info<3>: acc=0");
                auto const n4 = s.find("info<4>: acc=2");
                BOOST_TEST_NE(n1, s.npos);
                BOOST_TEST_NE(n2, s.npos);
                BOOST_TEST_NE(n3, s.npos);
                BOOST_TEST_NE(n4, s.npos);
                BOOST_TEST_EQ(counter, 4);
            }
            else
            BOOST_TEST_EQ(counter, 1);
#endif
        },
        []
        {
            std::abort();
        } );
    BOOST_TEST_EQ(counter, 0);

    leaf::try_handle_all(
        []() -> leaf::result<void>
        {
            return f2();
        },
        []
        {
            BOOST_TEST_EQ(counter, 0);
        } );

    return boost::report_errors();
}
