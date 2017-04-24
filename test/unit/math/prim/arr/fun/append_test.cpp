#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, append_array) {
  std::vector<double> x(3), y(2), z, result;

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  y[0] = 0.5;
  y[1] = 4.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0]);
  EXPECT_FLOAT_EQ(2.0, result[1]);
  EXPECT_FLOAT_EQ(3.0, result[2]);
  EXPECT_FLOAT_EQ(0.5, result[3]);
  EXPECT_FLOAT_EQ(4.0, result[4]);

  EXPECT_NO_THROW(result = stan::math::append_array(x, z));
  EXPECT_EQ(3, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0]);
  EXPECT_FLOAT_EQ(2.0, result[1]);
  EXPECT_FLOAT_EQ(3.0, result[2]);
}

TEST(MathFunctions, append_array_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> x(3), y(2), result;

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = nan;
  y[0] = 0.5;
  y[1] = 1.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_PRED1(boost::math::isnan<double>, result[2]);

  EXPECT_NO_THROW(result = stan::math::append_array(y, x));
  EXPECT_PRED1(boost::math::isnan<double>, result[4]);
}
