#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, isConstant) {
  using stan::is_constant;
  using stan::math::fvar;

  EXPECT_FALSE(is_constant<fvar<double> >::value);
}
