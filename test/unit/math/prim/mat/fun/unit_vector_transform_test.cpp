#include <stan/math/prim/mat.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

using Eigen::Dynamic;
using Eigen::Matrix;

TEST(prob_transform, unit_vector_rt0) {
  Matrix<double, Dynamic, 1> x(4);
  x << sqrt(0.1), -sqrt(0.2), -sqrt(0.3), sqrt(0.4);
  using stan::math::unit_vector_constrain;
  Matrix<double, Dynamic, 1> y = unit_vector_constrain(x);
  EXPECT_NEAR(x(0), y(0), 1e-8);
  EXPECT_NEAR(x(1), y(1), 1e-8);
  EXPECT_NEAR(x(2), y(2), 1e-8);
  EXPECT_NEAR(x(3), y(3), 1e-8);

  Matrix<double, Dynamic, 1> xrt = stan::math::unit_vector_free(y);
  EXPECT_EQ(x.size(), y.size());
  EXPECT_EQ(x.size(), xrt.size());
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_NEAR(x[i], xrt[i], 1E-10);
  }
}
TEST(prob_transform, unit_vector_rt) {
  Matrix<double, Dynamic, 1> x(3);
  x << 1.0, -1.0, 1.0;
  using stan::math::unit_vector_constrain;
  Matrix<double, Dynamic, 1> y = unit_vector_constrain(x);
  Matrix<double, Dynamic, 1> xrt = stan::math::unit_vector_free(y);
  EXPECT_EQ(x.size(), y.size());
  EXPECT_EQ(x.size(), xrt.size());
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_FLOAT_EQ(y[i], xrt[i]) << "error in component " << i;
  }
}
TEST(prob_transform, unit_vector_match) {
  Matrix<double, Dynamic, 1> x(3);
  x << 1.0, -1.0, 2.0;
  double lp;
  using stan::math::unit_vector_constrain;
  Matrix<double, Dynamic, 1> y = unit_vector_constrain(x);
  Matrix<double, Dynamic, 1> y2 = stan::math::unit_vector_constrain(x, lp);

  EXPECT_EQ(3, y.size());
  EXPECT_EQ(3, y2.size());
  for (int i = 0; i < x.size(); ++i)
    EXPECT_FLOAT_EQ(y[i], y2[i]) << "error in component " << i;
}
TEST(prob_transform, unit_vector_f_exception) {
  Matrix<double, Dynamic, 1> y(2);
  y << 0.5, 0.55;
  EXPECT_THROW(stan::math::unit_vector_free(y), std::domain_error);
  y << 1.1, -0.1;
  EXPECT_THROW(stan::math::unit_vector_free(y), std::domain_error);
  y << 0.0, 0.0;
  EXPECT_THROW(stan::math::unit_vector_free(y), std::domain_error);
  y(0) = std::numeric_limits<double>::quiet_NaN();
  y(1) = 1.0;
  EXPECT_THROW(stan::math::unit_vector_free(y), std::domain_error);
  y(0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::unit_vector_free(y), std::domain_error);
}
