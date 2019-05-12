#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, isConstantStruct) {
  using stan::is_constant_struct;

  EXPECT_FALSE(is_constant_struct<stan::math::var>::value);
  bool temp = is_constant_struct<stan::math::var, double>::value;
  EXPECT_FALSE(temp);
  temp = is_constant_struct<stan::math::var, double, stan::math::var>::value;
  EXPECT_FALSE(temp);
  temp = is_constant_struct<stan::math::var, double, double>::value;
  EXPECT_FALSE(temp);
  temp = is_constant_struct<stan::math::var, stan::math::var, stan::math::var, stan::math::var, stan::math::var,double, double>::value;
  EXPECT_FALSE(temp);
}
