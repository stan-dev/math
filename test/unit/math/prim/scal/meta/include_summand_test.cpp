#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::include_summand;

TEST(MetaTraits, IncludeSummandProptoFalse) {
  EXPECT_TRUE(include_summand<false>::value);
}

TEST(MetaTraits, IncludeSummandProptoTrueInt) {
  EXPECT_FALSE((include_summand<true, int>::value));
}

TEST(MetaTraits, IncludeSummandProptoTrueDouble) {
  EXPECT_FALSE((include_summand<true, double>::value));
}

TEST(MetaTraits, IncludeSummandConstantPropToTrueTen) {
  EXPECT_FALSE((include_summand<true, double, double, int, int, double, double,
                                int, int, double, int>::value));
}

TEST(MetaTraits, IncludeSummandConstantProptoFalseTen) {
  EXPECT_TRUE((include_summand<false, double, double, int, int, double, double,
                               int, int, double, int>::value));
}
