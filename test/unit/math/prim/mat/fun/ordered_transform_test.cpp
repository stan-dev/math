#include <stan/math/prim/mat.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using Eigen::Matrix;
using Eigen::Dynamic;

TEST(prob_transform,ordered) {
  Matrix<double,Dynamic,1> x(3);
  x << -15.0, -2.0, -5.0;
  Matrix<double,Dynamic,1> y = stan::math::ordered_constrain(x);
  EXPECT_EQ(x.size(), y.size());
  EXPECT_EQ(-15.0, y[0]);
  EXPECT_EQ(-15.0 + exp(-2.0), y[1]);
  EXPECT_EQ(-15.0 + exp(-2.0) + exp(-5.0), y[2]);
}
TEST(prob_transform,ordered_j) {
  Matrix<double,Dynamic,1> x(3);
  x << 1.0, -2.0, -5.0;
  double lp = -152.1;
  Matrix<double,Dynamic,1> y = stan::math::ordered_constrain(x,lp);
  EXPECT_EQ(x.size(), y.size());
  EXPECT_EQ(1.0, y[0]);
  EXPECT_EQ(1.0 + exp(-2.0), y[1]);
  EXPECT_EQ(1.0 + exp(-2.0) + exp(-5.0), y[2]);
  EXPECT_EQ(-152.1 - 2.0 - 5.0,lp);
}
TEST(prob_transform,ordered_f) {
  Matrix<double,Dynamic,1> y(3);
  y << -12.0, 1.1, 172.1;
  Matrix<double,Dynamic,1> x = stan::math::ordered_free(y);
  EXPECT_EQ(y.size(),x.size());
  EXPECT_FLOAT_EQ(-12.0, x[0]);
  EXPECT_FLOAT_EQ(log(1.1 + 12.0), x[1]);
  EXPECT_FLOAT_EQ(log(172.1 - 1.1), x[2]);
}
TEST(prob_transform,ordered_f_exception) {
  Matrix<double,Dynamic,1> y(3);
  y << -0.1, 0.0, 1.0;
  EXPECT_NO_THROW(stan::math::ordered_free(y));
  y << 0.0, 0.0, 0.0;
  EXPECT_THROW(stan::math::ordered_free(y), std::domain_error);
  y << 0.0, 1, 0.9;
  EXPECT_THROW(stan::math::ordered_free(y), std::domain_error);
}
TEST(prob_transform,ordered_rt) {
  Matrix<double,Dynamic,1> x(3);
  x << -1.0, 8.0, -3.9;
  Matrix<double,Dynamic,1> y = stan::math::ordered_constrain(x);
  Matrix<double,Dynamic,1> xrt = stan::math::ordered_free(y);
  EXPECT_EQ(x.size(), xrt.size());
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_FLOAT_EQ(x[i], xrt[i]);
  }
}
