#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::include_summand;
using stan::math::var;

TEST(MetaTraitsRevScal, IncludeSummandProptoTrueVar) {
  EXPECT_TRUE((include_summand<true, var>::value));
}

TEST(MetaTraitsRevScal, IncludeSummandProtoTrueVarTen) {
  EXPECT_TRUE((include_summand<true, double, var, int, var, double, double, int,
                               int, var, int>::value));
}
