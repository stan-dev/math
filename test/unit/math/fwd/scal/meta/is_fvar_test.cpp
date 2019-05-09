#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, is_fvar) {
  using stan::is_fvar;
  using stan::math::fvar;
  EXPECT_TRUE(is_fvar<fvar<int>>::value);
  EXPECT_TRUE(is_fvar<fvar<double>>::value);
  EXPECT_TRUE(is_fvar<fvar<fvar<double>>>::value);
  EXPECT_TRUE(is_fvar<fvar<fvar<fvar<double>>>>::value);
}
