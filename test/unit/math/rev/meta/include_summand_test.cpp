#include <stan/math/rev/meta.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using stan::math::include_summand;
using stan::math::var;

TEST(MetaTraitsRevScal, IncludeSummandProptoTrueVar) {
  EXPECT_TRUE((include_summand<true, var>::value));
}

TEST(MetaTraitsRevScal, IncludeSummandProtoTrueVarTen) {
  EXPECT_TRUE((include_summand<true, double, var, int, var, double, double, int,
                               int, var, int>::value));
}
