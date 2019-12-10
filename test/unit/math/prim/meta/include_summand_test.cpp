#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::include_summand;

TEST(MetaTraitsPrimScal, IncludeSummandProptoFalse) {
  EXPECT_TRUE(include_summand<false>::value);
}

TEST(MetaTraitsPrimScal, IncludeSummandProptoTrueInt) {
  EXPECT_FALSE((include_summand<true, int>::value));
}

TEST(MetaTraitsPrimScal, IncludeSummandProptoTrueDouble) {
  EXPECT_FALSE((include_summand<true, double>::value));
}

TEST(MetaTraitsPrimScal, IncludeSummandConstantPropToTrueTen) {
  EXPECT_FALSE((include_summand<true, double, double, int, int, double, double,
                                int, int, double, int>::value));
}

TEST(MetaTraitsPrimScal, IncludeSummandConstantProptoFalseTen) {
  EXPECT_TRUE((include_summand<false, double, double, int, int, double, double,
                               int, int, double, int>::value));
}
