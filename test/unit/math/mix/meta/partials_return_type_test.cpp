#include <stan/math/mix.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::partials_return_type;
using stan::math::fvar;
using stan::math::var;

TEST(MathMetaMix, PartialsReturnTypeFvarVar) {
  test::expect_same_type<var, partials_return_type<fvar<var> >::type>();
}

TEST(MathMetaMix, PartialsReturnTypeFvarFvarVar) {
  test::expect_same_type<fvar<var>,
                         partials_return_type<fvar<fvar<var> > >::type>();
}

TEST(MathMetaMix, PartialsReturnTypeFvarVarTenParams) {
  test::expect_same_type<
      var, partials_return_type<double, fvar<var>, double, int, double, float,
                                float, float, fvar<var>, int>::type>();
}

TEST(MathMetaMix, PartialsReturnTypeFvarFvarVarTenParams) {
  test::expect_same_type<
      fvar<var>,
      partials_return_type<double, fvar<fvar<var> >, double, int, double, float,
                           float, float, fvar<fvar<var> >, int>::type>();
}
