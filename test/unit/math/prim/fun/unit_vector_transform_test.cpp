#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, unit_vector_rt0) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  double lp = 0;
  Matrix<double, Dynamic, 1> x(4);
  x << sqrt(0.1), -sqrt(0.2), -sqrt(0.3), sqrt(0.4);
  std::vector<Matrix<double, Dynamic, 1>> x_vec{x, x, x};
  using stan::math::unit_vector_constrain;
  std::vector<Matrix<double, Dynamic, 1>> y_vec
      = unit_vector_constrain<false>(x_vec, lp);
  for (int i = 0; i < x_vec.size(); ++i) {
    EXPECT_EQ(x_vec[i].size(), y_vec[i].size());
    for (int j = 0; j < x_vec[i].size(); ++j) {
      EXPECT_NEAR(x_vec[i](j), y_vec[i](j), 1e-8);
    }
  }
  std::vector<Matrix<double, Dynamic, 1>> xrt
      = stan::math::unit_vector_free(y_vec);
  for (int i = 0; i < x_vec.size(); ++i) {
    EXPECT_EQ(x_vec[i].size(), xrt[i].size());
    for (int j = 0; j < x_vec[i].size(); ++j) {
      EXPECT_NEAR(x_vec[i](j), xrt[i](j), 1e-10);
    }
  }
  /*
  for (auto&& x_i : xrt) {
    for (int i = 0; i < x.size(); ++i) {
      EXPECT_NEAR(x[i], x_i[i], 1E-10);
    }
  }
  */
}
TEST(prob_transform, unit_vector_rt) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
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
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, 1> x(3);
  x << 1.0, -1.0, 2.0;
  double lp = 0;
  using stan::math::unit_vector_constrain;
  Matrix<double, Dynamic, 1> y = unit_vector_constrain(x);
  Matrix<double, Dynamic, 1> y2 = stan::math::unit_vector_constrain(x, lp);

  EXPECT_EQ(3, y.size());
  EXPECT_EQ(3, y2.size());
  for (int i = 0; i < x.size(); ++i)
    EXPECT_FLOAT_EQ(y[i], y2[i]) << "error in component " << i;
}
TEST(prob_transform, unit_vector_f_exception) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
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
