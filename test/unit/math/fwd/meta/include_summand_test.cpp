#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::fvar;
using stan::math::include_summand_b;

TEST(MathMetaFwd, IncludeSummandProptoTrueFvarDouble) {
  EXPECT_TRUE((include_summand_b<true, fvar<double>>));
}

TEST(MathMetaFwd, IncludeSummandProptoTrueFvarFvarDouble) {
  EXPECT_TRUE((include_summand_b<true, fvar<fvar<double>>>));
}

TEST(MathMetaFwd, IncludeSummandProtoTrueFvarDoubleTen) {
  EXPECT_TRUE((include_summand_b<true, double, fvar<double>, int, fvar<double>,
                                 double, double, int, int, fvar<double>, int>));
}

TEST(MathMetaFwd, IncludeSummandProtoTrueFvarFvarDoubleTen) {
  EXPECT_TRUE((include_summand_b<true, double, fvar<fvar<double>>, int,
                                 fvar<fvar<double>>, double, double, int, int,
                                 fvar<fvar<double>>, int>));
}
