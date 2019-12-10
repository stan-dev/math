#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, containsFvar) {
  using stan::contains_fvar;
  EXPECT_FALSE(contains_fvar<double>::value);
  bool temp
      = contains_fvar<double, double, double, double, double, double>::value;
  EXPECT_FALSE(temp);
  temp = contains_fvar<double, double, double, double, double, double, double,
                       double, double, double>::value;
  EXPECT_FALSE(temp);
}
