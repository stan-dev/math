#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(MathMetaPrim, IncludeSummandProptoFalse) {
  EXPECT_TRUE(stan::math::include_summand<false>::value);
}

TEST(MathMetaPrim, IncludeSummandProptoTrueInt) {
  EXPECT_FALSE((stan::math::include_summand<true, int>::value));
}

TEST(MathMetaPrim, IncludeSummandProptoTrueDouble) {
  EXPECT_FALSE((stan::math::include_summand<true, double>::value));
}

TEST(MathMetaPrim, IncludeSummandConstantPropToTrueTen) {
  EXPECT_FALSE(
      (stan::math::include_summand<true, double, double, int, int, double,
                                   double, int, int, double, int>::value));
}

TEST(MathMetaPrim, IncludeSummandConstantProptoFalseTen) {
  EXPECT_TRUE(
      (stan::math::include_summand<false, double, double, int, int, double,
                                   double, int, int, double, int>::value));
}
