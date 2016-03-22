#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, isConstantStruct) {
  using stan::is_constant_struct;

  EXPECT_FALSE(is_constant_struct<stan::math::var>::value);
}
