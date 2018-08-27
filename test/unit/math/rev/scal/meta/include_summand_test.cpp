#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::include_summand;
using stan::math::var;

TEST(MetaTraits, IncludeSummandProptoTrueVar) {
  EXPECT_TRUE((include_summand<true, var>::value));
}

TEST(MetaTraits, IncludeSummandProtoTrueVarTen) {
  EXPECT_TRUE((include_summand<true, double, var, int, var, double, double, int,
                               int, var, int>::value));
}
