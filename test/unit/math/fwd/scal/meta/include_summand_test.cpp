#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::fvar;
using stan::math::include_summand;

TEST(MetaTraits, IncludeSummandProptoTrueFvarDouble) {
  EXPECT_TRUE((include_summand<true, fvar<double> >::value));
}

TEST(MetaTraits, IncludeSummandProptoTrueFvarFvarDouble) {
  EXPECT_TRUE((include_summand<true, fvar<fvar<double> > >::value));
}

TEST(MetaTraits, IncludeSummandProtoTrueFvarDoubleTen) {
  EXPECT_TRUE(
      (include_summand<true, double, fvar<double>, int, fvar<double>, double,
                       double, int, int, fvar<double>, int>::value));
}

TEST(MetaTraits, IncludeSummandProtoTrueFvarFvarDoubleTen) {
  EXPECT_TRUE((include_summand<true, double, fvar<fvar<double> >, int,
                               fvar<fvar<double> >, double, double, int, int,
                               fvar<fvar<double> >, int>::value));
}
