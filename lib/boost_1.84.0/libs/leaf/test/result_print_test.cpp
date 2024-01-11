// Copyright 2018-2023 Emil Dotchevski and Reverge Studios, Inc.

// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef BOOST_LEAF_TEST_SINGLE_HEADER
#   include "leaf.hpp"
#else
#   include <boost/leaf/result.hpp>
#   include <boost/leaf/capture.hpp>
#   include <boost/leaf/context.hpp>
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
        std::stringstream ss;
        ss << r;
        std::string s = ss.str();
        std::cout << s << std::endl;
#if BOOST_LEAF_CFG_DIAGNOSTICS
        BOOST_TEST_EQ(s, "42");
#endif
#endif
    }

    {
        leaf::result<non_printable_value> r;
        BOOST_TEST(r);
#if BOOST_LEAF_CFG_STD_STRING
        std::stringstream ss;
        ss << r;
        std::string s = ss.str();
        std::cout << s << std::endl;
#if BOOST_LEAF_CFG_DIAGNOSTICS
        BOOST_TEST_EQ(s, "{not printable}");
#endif
#endif
    }

    {
        leaf::result<void> r;
        BOOST_TEST(r);
#if BOOST_LEAF_CFG_STD_STRING
        std::stringstream ss;
        ss << r;
        std::string s = ss.str();
        std::cout << s << std::endl;
#if BOOST_LEAF_CFG_DIAGNOSTICS
        BOOST_TEST_EQ(s, "No error");
#endif
#endif
    }

    {
        leaf::result<int> r = leaf::new_error(e_err{ });
        BOOST_TEST(!r);
#if BOOST_LEAF_CFG_STD_STRING
        std::stringstream ss;
        ss << r;
        std::string s = ss.str();
        std::cout << s << std::endl;
        leaf::error_id err = r.error();
#if BOOST_LEAF_CFG_DIAGNOSTICS
        BOOST_TEST_EQ(s, "Error ID " + std::to_string(err.value()));
#endif
#endif
    }

#if BOOST_LEAF_CFG_CAPTURE
    {
    using context_type = leaf::leaf_detail::polymorphic_context_impl<leaf::context<e_err>>;
        {
            leaf::result<int> r = leaf::capture( std::make_shared<context_type>(), []{ return leaf::result<int>( leaf::new_error( e_err { } ) ); } );
#if BOOST_LEAF_CFG_STD_STRING
            std::stringstream ss;
            ss << r;
            std::string s = ss.str();
            std::cout << s << std::endl;
            leaf::error_id err = r.error();
            BOOST_TEST_NE(s.find("Error ID " + std::to_string(err.value())), s.npos);
#if BOOST_LEAF_CFG_DIAGNOSTICS
            BOOST_TEST_NE(s.find("captured error objects"), s.npos);
            BOOST_TEST_NE(s.find("e_err"), s.npos);
#endif
#endif
        }
    }
#endif

    return boost::report_errors();
}
