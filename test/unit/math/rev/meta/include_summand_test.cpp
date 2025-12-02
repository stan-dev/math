#include <stan/math/rev/meta.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

using stan::math::include_summand;
using stan::math::var;

TEST_F(AgradRev, MetaTraitsRevScal_IncludeSummandProptoTrueVar) {
  EXPECT_TRUE((include_summand<true, var>::value));
}

TEST_F(AgradRev, MetaTraitsRevScal_IncludeSummandProtoTrueVarTen) {
  EXPECT_TRUE((include_summand<true, double, var, int, var, double, double, int,
                               int, var, int>::value));
}
