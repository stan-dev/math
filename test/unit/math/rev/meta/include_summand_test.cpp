#include <stan/math/rev/meta.hpp>
#include <test/unit/math/prim/scal/fun/promote_type_test_util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using stan::math::include_summand_b;
using stan::math::var;

TEST(MetaTraitsRevScal, IncludeSummandProptoTrueVar) {
  EXPECT_TRUE((include_summand_b<true, var>));
}

TEST(MetaTraitsRevScal, IncludeSummandProtoTrueVarTen) {
  EXPECT_TRUE((include_summand_b<true, double, var, int, var, double, double,
                                 int, int, var, int>));
}
