/*=============================================================================
    Copyright (c) 2001-2015 Joel de Guzman

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

    Duzy Chan:
    This test addresses the usage of boost::optional<foo> as an ast node.
=============================================================================*/
#include <boost/detail/lightweight_test.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/fusion/adapted/struct.hpp>
#include <boost/fusion/include/vector.hpp>

#include <iostream>
#include "test.hpp"
#include "utils.hpp"

struct twoints
{
    int a;
    int b;
};

struct adata
{
  twoints a;
  boost::optional<twoints> b;
  twoints c;
};

BOOST_FUSION_ADAPT_STRUCT(twoints, a, b)
BOOST_FUSION_ADAPT_STRUCT(adata, a, b, c)

int
main()
{
    {
      // Duzy Chan: This case addresses the usage of boost::optional<foo>
      // as an ast node. Which should actually test the ability of
      // boost::spirit::x3::traits::move_to to handle with optional source
      // value.
      boost::spirit::x3::rule<class twoints, adata> top = "top";
      boost::spirit::x3::rule<class twoints, boost::optional<twoints>>
        twoints = "twoints";
     
      using boost::spirit::x3::int_;
      auto const top_def = twoints >> ',' >> -twoints >> ',' >> twoints;
      auto const twoints_def = int_ >> int_;

      BOOST_SPIRIT_DEFINE(top, twoints);

      twoints a, b;
      BOOST_TEST((test_attr("1 2,3 4,5 6", top, a)));
      BOOST_TEST((a.a.a == 1 && a.a.b == 2));
      BOOST_TEST((a.b && a.b->a == 3 && a.b->b == 4));
      BOOST_TEST((a.c.a == 5 && a.c.b == 6));

      BOOST_TEST((test_attr("1 2,,5 6", top), b));
      BOOST_TEST((b.a.a == 1 && b.a.b == 2));
      BOOST_TEST((!a.b));
      BOOST_TEST((b.c.a == 5 && b.c.b == 6));
    }

    return boost::report_errors();
}
