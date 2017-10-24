#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>
#include <limits>

void test_unary_not(double x) {
  AVAR x_v = x;
  EXPECT_EQ(!x, !x_v);
}

TEST(AgradRev, unaryNot) {
  test_unary_not(6.1);
  test_unary_not(0);
  test_unary_not(-13.2);
  test_unary_not(std::numeric_limits<double>::infinity());
  test_unary_not(-std::numeric_limits<double>::infinity());
  test_unary_not(std::numeric_limits<double>::quiet_NaN());
}
