#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::fvar;
using stan::math::var;
using stan::return_type;

TEST(MathMetaMix, ReturnTypeFvarVar) {
  test::expect_same_type<fvar<var>, return_type<fvar<var> >::type>();
}

TEST(MathMetaMix, ReturnTypeFvarFvarVar) {
  test::expect_same_type<fvar<fvar<var> >,
                         return_type<fvar<fvar<var> > >::type>();
}

TEST(MathMetaMix, ReturnTypeFvarVarTenParams) {
  test::expect_same_type<
      fvar<var>, return_type<double, fvar<var>, double, int, double, float,
                             float, float, fvar<var>, int>::type>();
}

TEST(MathMetaMix, ReturnTypeFvarFvarVarTenParams) {
  test::expect_same_type<
      fvar<fvar<var> >,
      return_type<double, fvar<fvar<var> >, double, int, double, float, float,
                  float, fvar<fvar<var> >, int>::type>();
}
