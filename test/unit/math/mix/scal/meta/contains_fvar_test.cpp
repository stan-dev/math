#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, containsFvar) {
  using stan::contains_fvar;
  using stan::math::fvar;
  using stan::math::var;
  EXPECT_FALSE(contains_fvar<var>::value);
  EXPECT_TRUE((contains_fvar<double, fvar<var>, int>::value));
  bool temp
      = contains_fvar<double, fvar<var>, double, int, double, double>::value;
  EXPECT_TRUE(temp);
  temp = contains_fvar<double, double, int, double, double, double, int, double,
                       double, fvar<var>>::value;
  EXPECT_TRUE(temp);
}
