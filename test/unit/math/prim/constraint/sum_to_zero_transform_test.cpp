#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(prob_transform, sum_to_zero_rt0) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  double lp = 0;
  Matrix<double, Dynamic, 1> x(4);
  x << 0.0, 0.0, 0.0, 0.0;
  std::vector<Matrix<double, Dynamic, 1>> x_vec{x, x, x};
  std::vector<Matrix<double, Dynamic, 1>> y_vec
      = stan::math::sum_to_zero_constrain<false>(x_vec, lp);
  EXPECT_NO_THROW(stan::math::check_sum_to_zero("checkSumToZero", "y", y_vec));

  for (auto&& y_i : y_vec) {
    EXPECT_MATRIX_FLOAT_EQ(Eigen::VectorXd::Zero(5), y_i);
  }
  std::vector<Matrix<double, Dynamic, 1>> xrt
      = stan::math::sum_to_zero_free(y_vec);
  EXPECT_EQ(x.size() + 1, y_vec[2].size());
  for (auto&& x_i : xrt) {
    EXPECT_EQ(x.size(), x_i.size());
    for (int i = 0; i < x.size(); ++i) {
      EXPECT_NEAR(x[i], x_i[i], 1E-10);
    }
  }
}
TEST(prob_transform, sum_to_zero_rt) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, 1> x(3);
  x << 1.0, -1.0, 2.0;
  Matrix<double, Dynamic, 1> y = stan::math::sum_to_zero_constrain(x);
  EXPECT_NO_THROW(stan::math::check_sum_to_zero("checkSumToZero", "y", y));
  Matrix<double, Dynamic, 1> xrt = stan::math::sum_to_zero_free(y);
  EXPECT_EQ(x.size() + 1, y.size());
  EXPECT_EQ(x.size(), xrt.size());
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_FLOAT_EQ(x[i], xrt[i]);
  }
}
TEST(prob_transform, sum_to_zero_match) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, 1> x(3);
  x << 1.0, -1.0, 2.0;
  double lp = 0;
  Matrix<double, Dynamic, 1> y = stan::math::sum_to_zero_constrain(x);
  Matrix<double, Dynamic, 1> y2 = stan::math::sum_to_zero_constrain(x, lp);

  EXPECT_EQ(4, y.size());
  EXPECT_EQ(4, y2.size());
  for (int i = 0; i < y.size(); ++i)
    EXPECT_FLOAT_EQ(y[i], y2[i]);
}

TEST(prob_transform, sum_to_zero_f_exception) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, 1> y(2);
  y << 0.5, -0.55;
  EXPECT_THROW(stan::math::sum_to_zero_free(y), std::domain_error);
}

TEST(prob_transform, sum_to_zero_mat_rt0) {
  double lp = 0;
  auto x = Eigen::MatrixXd::Zero(4, 3);
  std::vector<Eigen::MatrixXd> x_vec{x, x, x};
  std::vector<Eigen::MatrixXd> y_vec
      = stan::math::sum_to_zero_constrain<false>(x_vec, lp);
  EXPECT_NO_THROW(stan::math::check_sum_to_zero("checkSumToZero", "y", y_vec));

  for (auto&& y_i : y_vec) {
    EXPECT_MATRIX_FLOAT_EQ(Eigen::MatrixXd::Zero(5, 4), y_i);
  }
  std::vector<Eigen::MatrixXd> xrt = stan::math::sum_to_zero_free(y_vec);
  EXPECT_EQ(x.rows() + 1, y_vec[2].rows());
  EXPECT_EQ(x.cols() + 1, y_vec[2].cols());
  for (auto&& x_i : xrt) {
    EXPECT_MATRIX_NEAR(x, x_i, 1E-10);
  }
}

TEST(prob_transform, sum_to_zero_mat_rt) {
  Eigen::MatrixXd x(3, 4);
  x << 1.0, -1.0, 2.0, 1.1, -1.1, 2.1, 1.2, -1.2, 2.2, 1.3, -1.3, 2.3;

  Eigen::MatrixXd y = stan::math::sum_to_zero_constrain(x);

  EXPECT_EQ(x.rows() + 1, y.rows());
  EXPECT_EQ(x.cols() + 1, y.cols());
  EXPECT_FLOAT_EQ(y.sum(), 0.0);
  EXPECT_MATRIX_NEAR(y.rowwise().sum(), Eigen::VectorXd::Zero(4), 1E-15);
  EXPECT_MATRIX_NEAR(y.colwise().sum(), Eigen::RowVectorXd::Zero(5), 1E-15);
  EXPECT_NO_THROW(stan::math::check_sum_to_zero("checkSumToZero", "y", y));

  Eigen::MatrixXd xrt = stan::math::sum_to_zero_free(y);

  EXPECT_EQ(x.size(), xrt.size());
  EXPECT_MATRIX_FLOAT_EQ(x, xrt);
}

TEST(prob_transform, sum_to_zero_mat_match) {
  Eigen::MatrixXd x(3, 4);
  x << 1.0, -1.0, 2.0, 1.1, -1.1, 2.1, 1.2, -1.2, 2.2, 1.3, -1.3, 2.3;

  double lp = 0;
  Eigen::MatrixXd y = stan::math::sum_to_zero_constrain(x);
  Eigen::MatrixXd y2 = stan::math::sum_to_zero_constrain(x, lp);

  EXPECT_EQ(y.size(), y2.size());
  EXPECT_MATRIX_FLOAT_EQ(y, y2);
  EXPECT_EQ(lp, 0);
}

TEST(prob_transform, sum_to_zero_mat_f_exception) {
  Eigen::MatrixXd x(3, 4);
  x << 1.0, -1.0, 2.0, 1.1, -1.1, 2.1, 1.2, -1.2, 2.2, 1.3, -1.3, 2.3;
  auto y = stan::math::sum_to_zero_constrain(x);
  y(0, 0) += 0.1;
  EXPECT_THROW(stan::math::sum_to_zero_free(y), std::domain_error);
}
