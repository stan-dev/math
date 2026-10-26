#include <stan/math/mix.hpp>
#include <gtest/gtest.h>

TEST(MathMetaMix, IncludeSummandProptoTrueFvarVar) {
  using stan::math::fvar;
  using stan::math::include_summand;
  using stan::math::var;
  EXPECT_TRUE((include_summand<true, fvar<var> >::value));
  EXPECT_TRUE((include_summand<true, fvar<fvar<var> > >::value));
  EXPECT_TRUE((include_summand<true, double, fvar<var>, int, fvar<var>, double,
                               double, int, int, fvar<var>, int>::value));
  EXPECT_TRUE((
      include_summand<true, double, fvar<fvar<var> >, int, fvar<fvar<var> >,
                      double, double, int, int, fvar<fvar<var> >, int>::value));
}
