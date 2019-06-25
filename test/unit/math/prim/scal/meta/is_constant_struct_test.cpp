#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, isConstantStruct) {
  using stan::is_constant_struct;
  EXPECT_TRUE(is_constant_all<int>::value);
  EXPECT_TRUE(is_constant_all<double>::value);
  EXPECT_TRUE(is_constant_all<float>::value);
  EXPECT_TRUE(is_constant_all<int32_t>::value);
  bool temp = is_constant_all<int, int>::value;
  EXPECT_TRUE(temp);
  temp = is_constant_all<double, double, double>::value;
  EXPECT_TRUE(temp);
  temp = is_constant_all<float, float, float, float>::value;
  EXPECT_TRUE(temp);
  temp = is_constant_all<int32_t, int32_t, int32_t, int32_t>::value;
  EXPECT_TRUE(temp);
}
