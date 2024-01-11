// Copyright 2018-2023 Emil Dotchevski and Reverge Studios, Inc.

// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/leaf/config.hpp>

#if defined(BOOST_LEAF_NO_EXCEPTIONS) || !BOOST_LEAF_CFG_CAPTURE

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
#   include <boost/leaf/capture.hpp>
#   include <boost/leaf/handle_errors.hpp>
#   include <boost/leaf/exception.hpp>
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
    auto error_handlers = std::make_tuple(
        []( info<1>, info<3> )
        {
            return 1;
        },
        []
        {
            return 2;
        } );
    BOOST_TEST_EQ(count, 0);
    std::exception_ptr ep;
    leaf::context_ptr ctx = leaf::make_shared_context(error_handlers);
    try
    {
        leaf::capture(
            ctx,
            []
            {
                leaf::throw_exception(info<1>{}, info<3>{});
            } );
        BOOST_TEST(false);
    }
    catch(...)
    {
        ep = std::current_exception();
    }
    BOOST_TEST_EQ(count, 2);

#if BOOST_LEAF_CFG_STD_STRING
    {
        std::ostringstream st;
        st << *ctx;
        std::string s = st.str();
        std::cout << s << std::endl;
#if BOOST_LEAF_CFG_DIAGNOSTICS
        BOOST_TEST_NE(s.find("info<1> instance"), s.npos);
        BOOST_TEST_NE(s.find("info<3> instance"), s.npos);
#endif
    }
#endif

    int r = leaf::try_catch(
        [&]() -> int
        {
            std::rethrow_exception(ep);
        },
        error_handlers );
    BOOST_TEST_EQ(r, 1);
    ep = std::exception_ptr();
    BOOST_TEST_EQ(count, 0);

    return boost::report_errors();
}

#endif
