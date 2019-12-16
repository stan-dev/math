#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::fvar;
using stan::math::include_summand;
using stan::math::include_summand_b;
using stan::math::var;

TEST(MathMetaMix, IncludeSummandProptoTrueFvarVar) {
  EXPECT_TRUE((include_summand<true, fvar<var> >::value));
  EXPECT_TRUE((include_summand_b<true, fvar<var>>));
}

TEST(MathMetaMix, IncludeSummandProptoTrueFvarFvarVar) {
  EXPECT_TRUE((include_summand<true, fvar<fvar<var> > >::value));
  EXPECT_TRUE((include_summand_b<true, fvar<fvar<var>>>));
}

TEST(MathMetaMix, IncludeSummandProtoTrueFvarVarTen) {
  EXPECT_TRUE((include_summand<true, double, fvar<var>, int, fvar<var>, double,
                               double, int, int, fvar<var>, int>::value));
  EXPECT_TRUE((include_summand_b<true, double, fvar<var>, int, fvar<var>,
                                 double, double, int, int, fvar<var>, int>));
}

TEST(MathMetaMix, IncludeSummandProtoTrueFvarFvarVarTen) {
  EXPECT_TRUE((
      include_summand<true, double, fvar<fvar<var> >, int, fvar<fvar<var> >,
                      double, double, int, int, fvar<fvar<var> >, int>::value));
  EXPECT_TRUE((
      include_summand_b<true, double, fvar<fvar<var>>, int, fvar<fvar<var>>,
                        double, double, int, int, fvar<fvar<var>>, int>));
}
