
#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>




TEST(MetaTraits_scal, isConstant) {
  using stan::is_constant;

  EXPECT_TRUE(is_constant<double>::value);
  EXPECT_TRUE(is_constant<float>::value);
  EXPECT_TRUE(is_constant<unsigned int>::value);
  EXPECT_TRUE(is_constant<int>::value);
}
