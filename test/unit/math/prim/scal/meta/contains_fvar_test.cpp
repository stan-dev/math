#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, containsFvar) {
  using stan::contains_fvar;
  EXPECT_FALSE(contains_fvar<double>::value);
  bool temp
      = contains_fvar<double, double, double, double, double, double>::value;
  EXPECT_FALSE(temp);
  temp = contains_fvar<double, double, double, double, double, double, double,
                       double, double, double>::value;
  EXPECT_FALSE(temp);
}
