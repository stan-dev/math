#include <stan/math/prim/arr.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/fwd/core/fvar.hpp>
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

  EXPECT_NO_THROW(result = stan::math::append_array(z, z));
  EXPECT_EQ(0, result.size());
}

TEST(MathFunctions, append_array_int) {
  std::vector<double> x(3), result;
  std::vector<int> y(2), z(3);

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  y[0] = 5;
  y[1] = 4;
  z[0] = 6;
  z[1] = 7;
  z[2] = 8;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0]);
  EXPECT_FLOAT_EQ(2.0, result[1]);
  EXPECT_FLOAT_EQ(3.0, result[2]);
  EXPECT_FLOAT_EQ(5.0, result[3]);
  EXPECT_FLOAT_EQ(4.0, result[4]);

  EXPECT_NO_THROW(result = stan::math::append_array(y, z));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(5.0, result[0]);
  EXPECT_FLOAT_EQ(4.0, result[1]);
  EXPECT_FLOAT_EQ(6.0, result[2]);
  EXPECT_FLOAT_EQ(7.0, result[3]);
  EXPECT_FLOAT_EQ(8.0, result[4]);
}

TEST(MathFunctions, append_array_var) {
  std::vector<double> x(3);
  std::vector<stan::math::var> y(2), result;

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  y[0] = 0.5;
  y[1] = 4.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0].val());
  EXPECT_FLOAT_EQ(2.0, result[1].val());
  EXPECT_FLOAT_EQ(3.0, result[2].val());
  EXPECT_FLOAT_EQ(0.5, result[3].val());
  EXPECT_FLOAT_EQ(4.0, result[4].val());
}

TEST(MathFunctions, append_array_fvar) {
  std::vector<double> x(3);
  std::vector<stan::math::fvar<double> > y(2), result;

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  y[0] = 0.5;
  y[1] = 4.0;

  EXPECT_NO_THROW(result = stan::math::append_array(y, x));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(0.5, result[0].val());
  EXPECT_FLOAT_EQ(4.0, result[1].val());
  EXPECT_FLOAT_EQ(1.0, result[2].val());
  EXPECT_FLOAT_EQ(2.0, result[3].val());
  EXPECT_FLOAT_EQ(3.0, result[4].val());
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
