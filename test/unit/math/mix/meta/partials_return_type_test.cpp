#include <stan/math/mix.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(MathMetaMix, PartialsReturnTypeFvarVar) {
  using stan::partials_return_type;
  using stan::math::fvar;
  using stan::math::var;
  EXPECT_SAME_TYPE(var, partials_return_type<fvar<var> >::type);
}

TEST(MathMetaMix, PartialsReturnTypeFvarFvarVar) {
  using stan::partials_return_type;
  using stan::math::fvar;
  using stan::math::var;
  EXPECT_SAME_TYPE(fvar<var>, partials_return_type<fvar<fvar<var> > >::type);
}

TEST(MathMetaMix, PartialsReturnTypeFvarVarTenParams) {
  using stan::partials_return_type;
  using stan::math::fvar;
  using stan::math::var;
  EXPECT_SAME_TYPE(
      var, partials_return_type<double, fvar<var>, double, int, double, float,
                                float, float, fvar<var>, int>::type);
}

TEST(MathMetaMix, PartialsReturnTypeFvarFvarVarTenParams) {
  using stan::partials_return_type;
  using stan::math::fvar;
  using stan::math::var;
  EXPECT_SAME_TYPE(
      fvar<var>,
      partials_return_type<double, fvar<fvar<var> >, double, int, double, float,
                           float, float, fvar<fvar<var> >, int>::type);
}
