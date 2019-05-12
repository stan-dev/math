#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, isConstantStruct) {
  using stan::is_constant_struct;
  EXPECT_TRUE(is_constant_struct<int>::value);
  EXPECT_TRUE(is_constant_struct<double>::value);
  EXPECT_TRUE(is_constant_struct<float>::value);
  EXPECT_TRUE(is_constant_struct<int32_t>::value);
  bool temp = is_constant_struct<int, int>::value;
  EXPECT_TRUE(temp);
  temp = is_constant_struct<double, double, double>::value;
  EXPECT_TRUE(temp);
  temp = is_constant_struct<float, float, float, float>::value;
  EXPECT_TRUE(temp);
  temp = is_constant_struct<int32_t, int32_t, int32_t, int32_t>::value;
  EXPECT_TRUE(temp);
}
