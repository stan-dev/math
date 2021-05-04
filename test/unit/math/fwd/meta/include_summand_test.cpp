#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::fvar;
using stan::math::include_summand;

TEST(MathMetaFwd, IncludeSummandProptoTrueFvarDouble) {
  EXPECT_TRUE((include_summand<true, fvar<double> >::value));
}

TEST(MathMetaFwd, IncludeSummandProptoTrueFvarFvarDouble) {
  EXPECT_TRUE((include_summand<true, fvar<fvar<double> > >::value));
}

TEST(MathMetaFwd, IncludeSummandProtoTrueFvarDoubleTen) {
  EXPECT_TRUE(
      (include_summand<true, double, fvar<double>, int, fvar<double>, double,
                       double, int, int, fvar<double>, int>::value));
}

TEST(MathMetaFwd, IncludeSummandProtoTrueFvarFvarDoubleTen) {
  EXPECT_TRUE((include_summand<true, double, fvar<fvar<double> >, int,
                               fvar<fvar<double> >, double, double, int, int,
                               fvar<fvar<double> >, int>::value));
}
