#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(prob_transform, simplex_rt0) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  double lp = 0;
  Matrix<double, Dynamic, 1> x(4);
  x << 0.0, 0.0, 0.0, 0.0;
  std::vector<Matrix<double, Dynamic, 1>> x_vec{x, x, x};
  std::vector<Matrix<double, Dynamic, 1>> y_vec
      = stan::math::simplex_constrain<false>(x_vec, lp);
  for (auto&& y_i : y_vec) {
    EXPECT_MATRIX_FLOAT_EQ(Eigen::VectorXd::Constant(5, 1.0 / 5.0), y_i);
  }
  std::vector<Matrix<double, Dynamic, 1>> xrt = stan::math::simplex_free(y_vec);
  EXPECT_EQ(x.size() + 1, y_vec[2].size());
  for (auto&& x_i : xrt) {
    EXPECT_EQ(x.size(), x_i.size());
    for (int i = 0; i < x.size(); ++i) {
      EXPECT_NEAR(x[i], x_i[i], 1E-10);
    }
  }
}
TEST(prob_transform, simplex_rt) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
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
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, 1> x(3);
  x << 1.0, -1.0, 2.0;
  double lp = 0;
  Matrix<double, Dynamic, 1> y = stan::math::simplex_constrain(x);
  Matrix<double, Dynamic, 1> y2 = stan::math::simplex_constrain(x, lp);

  EXPECT_EQ(4, y.size());
  EXPECT_EQ(4, y2.size());
  for (int i = 0; i < x.size(); ++i)
    EXPECT_FLOAT_EQ(y[i], y2[i]);
}
TEST(prob_transform, simplex_f_exception) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, 1> y(2);
  y << 0.5, 0.55;
  EXPECT_THROW(stan::math::simplex_free(y), std::domain_error);
  y << 1.1, -0.1;
  EXPECT_THROW(stan::math::simplex_free(y), std::domain_error);
}
