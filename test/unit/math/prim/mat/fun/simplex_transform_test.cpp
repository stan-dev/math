#include <stan/math/prim/mat.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using Eigen::Dynamic;
using Eigen::Matrix;

TEST(prob_transform, simplex_rt0) {
  Matrix<double, Dynamic, 1> x(4);
  x << 0.0, 0.0, 0.0, 0.0;
  Matrix<double, Dynamic, 1> y = stan::math::simplex_constrain(x);
  EXPECT_FLOAT_EQ(1.0 / 5.0, y(0));
  EXPECT_FLOAT_EQ(1.0 / 5.0, y(1));
  EXPECT_FLOAT_EQ(1.0 / 5.0, y(2));
  EXPECT_FLOAT_EQ(1.0 / 5.0, y(3));
  EXPECT_FLOAT_EQ(1.0 / 5.0, y(4));

  Matrix<double, Dynamic, 1> xrt = stan::math::simplex_free(y);
  EXPECT_EQ(x.size() + 1, y.size());
  EXPECT_EQ(x.size(), xrt.size());
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_NEAR(x[i], xrt[i], 1E-10);
  }
}
TEST(prob_transform, simplex_rt) {
  Matrix<double, Dynamic, 1> x(3);
  x << 1.0, -1.0, 2.0;
  Matrix<double, Dynamic, 1> y = stan::math::simplex_constrain(x);
  Matrix<double, Dynamic, 1> xrt = stan::math::simplex_free(y);
  EXPECT_EQ(x.size() + 1, y.size());
  EXPECT_EQ(x.size(), xrt.size());
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_FLOAT_EQ(x[i], xrt[i]);
  }
}
TEST(prob_transform, simplex_match) {
  Matrix<double, Dynamic, 1> x(3);
  x << 1.0, -1.0, 2.0;
  double lp;
  Matrix<double, Dynamic, 1> y = stan::math::simplex_constrain(x);
  Matrix<double, Dynamic, 1> y2 = stan::math::simplex_constrain(x, lp);

  EXPECT_EQ(4, y.size());
  EXPECT_EQ(4, y2.size());
  for (int i = 0; i < x.size(); ++i)
    EXPECT_FLOAT_EQ(y[i], y2[i]);
}
TEST(prob_transform, simplex_f_exception) {
  Matrix<double, Dynamic, 1> y(2);
  y << 0.5, 0.55;
  EXPECT_THROW(stan::math::simplex_free(y), std::domain_error);
  y << 1.1, -0.1;
  EXPECT_THROW(stan::math::simplex_free(y), std::domain_error);
}
