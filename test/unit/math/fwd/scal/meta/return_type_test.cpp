#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::fvar;
using stan::return_type;

TEST(MetaTraits, ReturnTypeFvar) {
  test::expect_same_type<fvar<double>, return_type<fvar<double> >::type>();
  test::expect_same_type<fvar<double>,
                         return_type<fvar<double>, fvar<double> >::type>();
  test::expect_same_type<fvar<double>, return_type<fvar<double>, fvar<double>,
                                                   fvar<double> >::type>();
  test::expect_same_type<
      fvar<double>, return_type<fvar<double>, fvar<double>, double>::type>();
  test::expect_same_type<fvar<fvar<double> >,
                         return_type<fvar<fvar<double> > >::type>();
  test::expect_same_type<
      fvar<fvar<double> >,
      return_type<fvar<fvar<double> >, fvar<fvar<double> > >::type>();
  test::expect_same_type<fvar<fvar<double> >,
                         return_type<fvar<fvar<double> >, fvar<fvar<double> >,
                                     fvar<fvar<double> > >::type>();
  test::expect_same_type<
      fvar<fvar<double> >,
      return_type<fvar<fvar<double> >, fvar<fvar<double> >, double>::type>();
}
