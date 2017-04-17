#include <stan/math/prim/scal.hpp>
#include <test/unit/math/prim/scal/fun/promote_type_test_util.hpp>
#include <gtest/gtest.h>

TEST(MathFunctionsPromoteScalar, Match) {
  using stan::math::promote_scalar;
  EXPECT_FLOAT_EQ(1.3, promote_scalar<double>(1.3));
  EXPECT_EQ(3, promote_scalar<int>(3));

  expect_type<double>(promote_scalar<double>(2.3));
  expect_type<int>(promote_scalar<int>(2));
}
TEST(MathFunctionsPromoteScalar, Mismatch) {
  using stan::math::promote_scalar;
  EXPECT_FLOAT_EQ(2.0, promote_scalar<double>(2));
  expect_type<double>(promote_scalar<double>(2));
}
