
#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>





using stan::math::include_summand;

TEST(MetaTraits_scal, IncludeSummandProptoFalse) {
  EXPECT_TRUE(include_summand<false>::value);
}

TEST(MetaTraits_scal, IncludeSummandProptoTrueInt) {
  EXPECT_FALSE((include_summand<true, int>::value));
}

TEST(MetaTraits_scal, IncludeSummandProptoTrueDouble) {
  EXPECT_FALSE((include_summand<true, double>::value));
}

TEST(MetaTraits_scal, IncludeSummandConstantPropToTrueTen) {
  EXPECT_FALSE((include_summand<true, double, double, int, int, double, double,
                                int, int, double, int>::value));
}

TEST(MetaTraits_scal, IncludeSummandConstantProptoFalseTen) {
  EXPECT_TRUE((include_summand<false, double, double, int, int, double, double,
                               int, int, double, int>::value));
}
