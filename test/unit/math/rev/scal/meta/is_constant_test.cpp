#include <gtest/gtest.h>
#include <stan/math/rev/scal.hpp>

TEST(MetaTraits, isConstant) {
  using stan::is_constant;
  using stan::math::var;

  EXPECT_FALSE(is_constant<var>::value);
}
