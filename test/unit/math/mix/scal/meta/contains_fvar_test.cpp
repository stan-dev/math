#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits,containsFvar) {
  using stan::math::var;
  using stan::math::fvar;
  using stan::contains_fvar;
  EXPECT_FALSE(contains_fvar<var>::value);
  EXPECT_TRUE((contains_fvar<double, fvar<var>, int >::value));
}
