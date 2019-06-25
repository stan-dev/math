#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, isConstant) {
  using stan::is_constant;

  EXPECT_TRUE(is_constant_all<double>::value);
  EXPECT_TRUE(is_constant_all<float>::value);
  EXPECT_TRUE(is_constant_all<unsigned int>::value);
  EXPECT_TRUE(is_constant_all<int>::value);
}
