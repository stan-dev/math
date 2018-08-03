#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

void test_unary_not(double x) {
  using stan::math::fvar;
  using stan::math::var;

  EXPECT_EQ(!x, !fvar<double>(x));
  EXPECT_EQ(!x, !fvar<fvar<double> >(x));
  EXPECT_EQ(!x, !fvar<var>(x));
  EXPECT_EQ(!x, !fvar<fvar<var> >(x));
}

TEST(AgradRev, unaryNot) {
  test_unary_not(6.1);
  test_unary_not(0);
  test_unary_not(-13.2);
  test_unary_not(std::numeric_limits<double>::infinity());
  test_unary_not(-std::numeric_limits<double>::infinity());
  test_unary_not(std::numeric_limits<double>::quiet_NaN());
}
