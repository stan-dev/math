#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::fvar;
using stan::return_type;

TEST(MathMetaFwd, ReturnTypeFvarDouble) {
  test::expect_same_type<fvar<double>, return_type<fvar<double> >::type>();
}

TEST(MathMetaFwd, ReturnTypeFvarFvarDouble) {
  test::expect_same_type<fvar<fvar<double> >,
                         return_type<fvar<fvar<double> > >::type>();
}

TEST(MathMetaFwd, ReturnTypeFvarDoubleTenParams) {
  test::expect_same_type<
      fvar<double>,
      return_type<double, fvar<double>, double, int, double, float, float,
                  float, fvar<double>, int>::type>();
}

TEST(MathMetaFwd, ReturnTypeFvarFvarDoubleTenParams) {
  test::expect_same_type<
      fvar<fvar<double> >,
      return_type<double, fvar<fvar<double> >, double, int, double, float,
                  float, float, fvar<fvar<double> >, int>::type>();
}
