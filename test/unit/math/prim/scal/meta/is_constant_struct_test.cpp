#include <gtest/gtest.h>
#include <stan/math/prim/scal.hpp>

TEST(MetaTraits, isConstantStruct) {
  using stan::is_constant_struct;
  EXPECT_TRUE(is_constant_struct<int>::value);
  EXPECT_TRUE(is_constant_struct<double>::value);
  EXPECT_TRUE(is_constant_struct<float>::value);
  EXPECT_TRUE(is_constant_struct<int32_t>::value);
}
