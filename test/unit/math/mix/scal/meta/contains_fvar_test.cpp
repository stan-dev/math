#include <gtest/gtest.h>
#include <stan/math/mix/scal.hpp>

TEST(MetaTraits, containsFvar) {
  using stan::contains_fvar;
  using stan::math::fvar;
  using stan::math::var;
  EXPECT_FALSE(contains_fvar<var>::value);
  EXPECT_TRUE((contains_fvar<double, fvar<var>, int>::value));
}
