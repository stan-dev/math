// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef BOOST_LEAF_TEST_SINGLE_HEADER
#   include "leaf.hpp"
#else
#   include <boost/leaf/detail/demangle.hpp>
#endif

#include <cstring>
#include <cstdint>

namespace leaf = boost::leaf;

#include "lightweight_test.hpp"

namespace leaf_test
{
    class class_ { };
    struct struct_ { };
    enum enum_ { };
    template <int> class class_template1 { };
    template <int> struct struct_template1 { };
    template <class> class class_template2 { };
    template <class> struct struct_template2 { };
}

namespace boost { namespace leaf {
    struct in_namespace_boost_leaf { };
} }

class class_ { };
struct struct_ { };
enum enum_ { };
template <int> class class_template1 { };
template <int> struct struct_template1 { };
template <class> class class_template2 { };
template <class> struct struct_template2 { };

bool test(leaf::parsed const & pn, char const * correct)
{
    return
        std::strlen(correct) == pn.len &&
        std::memcmp(correct, pn.name, pn.len) == 0;
}

int main()
{
    using leaf::parse;

    BOOST_TEST(test(parse<leaf::in_namespace_boost_leaf>(), "boost::leaf::in_namespace_boost_leaf"));

    BOOST_TEST(test(parse<int>(), "int"));

    BOOST_TEST(test(parse<leaf_test::class_>(), "leaf_test::class_"));
    BOOST_TEST(test(parse<leaf_test::struct_>(), "leaf_test::struct_"));
    BOOST_TEST(test(parse<leaf_test::enum_>(), "leaf_test::enum_"));
    BOOST_TEST(test(parse<leaf_test::class_template1<42>>(), "leaf_test::class_template1<42>"));
    BOOST_TEST(test(parse<leaf_test::struct_template1<42>>(), "leaf_test::struct_template1<42>"));
    BOOST_TEST(test(parse<leaf_test::class_template2<int>>(), "leaf_test::class_template2<int>"));
    BOOST_TEST(test(parse<leaf_test::struct_template2<int>>(), "leaf_test::struct_template2<int>"));

    BOOST_TEST(test(parse<class_>(), "class_"));
    BOOST_TEST(test(parse<struct_>(), "struct_"));
    BOOST_TEST(test(parse<enum_>(), "enum_"));
    BOOST_TEST(test(parse<class_template1<42>>(), "class_template1<42>"));
    BOOST_TEST(test(parse<struct_template1<42>>(), "struct_template1<42>"));
    BOOST_TEST(test(parse<class_template2<int>>(), "class_template2<int>"));
    BOOST_TEST(test(parse<struct_template2<int>>(), "struct_template2<int>"));

    return boost::report_errors();
}
