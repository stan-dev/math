#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::include_summand_b;

TEST(MathMetaPrim, IncludeSummandProptoFalse) {
  EXPECT_TRUE(include_summand_b<false>);
}

TEST(MathMetaPrim, IncludeSummandProptoTrueInt) {
  EXPECT_FALSE((include_summand_b<true, int>));
}

TEST(MathMetaPrim, IncludeSummandProptoTrueDouble) {
  EXPECT_FALSE((include_summand_b<true, double>));
}

TEST(MathMetaPrim, IncludeSummandConstantPropToTrueTen) {
  EXPECT_FALSE((include_summand_b<true, double, double, int, int, double,
                                  double, int, int, double, int>));
}

TEST(MathMetaPrim, IncludeSummandConstantProptoFalseTen) {
  EXPECT_TRUE((include_summand_b<false, double, double, int, int, double,
                                 double, int, int, double, int>));
}
