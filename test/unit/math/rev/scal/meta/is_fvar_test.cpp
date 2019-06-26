#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, is_fvar) {
  using stan::is_fvar;
  EXPECT_FALSE(is_fvar<stan::math::var>::value);
}
