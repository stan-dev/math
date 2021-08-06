#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, checkCorrCholeskyMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_mat;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y_vec{
      y_mat, y_mat, y_mat};

  using stan::math::check_cholesky_factor_corr;
  using std::sqrt;
  for (auto& y : y_vec) {
    y.resize(1, 1);
    y << 1;
  }
  EXPECT_NO_THROW(
      check_cholesky_factor_corr("checkCorrCholeskyMatrix", "y", y_vec));

  for (auto& y : y_vec) {
    y.resize(3, 3);
    y << 1, 0, 0, sqrt(0.5), sqrt(0.5), 0, sqrt(0.25), sqrt(0.25), sqrt(0.5);
  }
  EXPECT_NO_THROW(
      check_cholesky_factor_corr("checkCorrCholeskyMatrix", "y", y_vec));

  // not positive
  for (auto& y : y_vec) {
    y.resize(1, 1);
    y << -1;
  }
  EXPECT_THROW(
      check_cholesky_factor_corr("checkCorrCholeskyMatrix", "y", y_vec),
      std::domain_error);

  // not lower triangular
  for (auto& y : y_vec) {
    y.resize(3, 3);
    y << 1, 2, 3, 0, 5, 6, 0, 0, 9;
  }
  EXPECT_THROW(
      check_cholesky_factor_corr("checkCorrCholeskyMatrix", "y", y_vec),
      std::domain_error);

  // not positive
  for (auto& y : y_vec) {
    y.resize(3, 3);
    y << 1, 0, 0, 2, -1, 0, 1, 2, 3;
  }
  EXPECT_THROW(
      check_cholesky_factor_corr("checkCorrCholeskyMatrix", "y", y_vec),
      std::domain_error);

  // not rectangular
  for (auto& y : y_vec) {
    y.resize(2, 3);
    y << 1, 2, 3, 4, 5, 6;
  }
  EXPECT_THROW(
      check_cholesky_factor_corr("checkCorrCholeskyMatrix", "y", y_vec),
      std::invalid_argument);
  for (auto& y : y_vec) {
    y.resize(3, 2);
    y << 1, 0, 2, 3, 4, 5;
  }
  EXPECT_THROW(
      check_cholesky_factor_corr("checkCorrCholeskyMatrix", "y", y_vec),
      std::invalid_argument);

  // not unit vectors
  for (auto& y : y_vec) {
    y.resize(3, 3);
    y << 1, 0, 0, 1, 1, 0, 1, 1, 1;
  }
  EXPECT_THROW(
      check_cholesky_factor_corr("checkCorrCholeskyMatrix", "y", y_vec),
      std::domain_error);
}

TEST(ErrorHandlingMatrix, checkCorrCholeskyMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_mat;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y_vec{
      y_mat, y_mat, y_mat};
  double nan = std::numeric_limits<double>::quiet_NaN();

  using stan::math::check_cholesky_factor_corr;
  using std::sqrt;

  for (auto& y : y_vec) {
    y.resize(1, 1);
    y << nan;
  }
  EXPECT_THROW(
      check_cholesky_factor_corr("checkCorrCholeskyMatrix", "y", y_vec),
      std::domain_error);

  for (auto& y : y_vec) {
    y.resize(3, 3);
    y << 1, 0, 0, sqrt(0.5), sqrt(0.5), 0, sqrt(0.25), sqrt(0.25), sqrt(0.5);
  }
  EXPECT_NO_THROW(
      check_cholesky_factor_corr("checkCorrCholeskyMatrix", "y", y_vec));

  for (int i = 0; i < y_vec[0].size(); i++) {
    for (auto& y : y_vec) {
      y(i) = nan;
    }
    EXPECT_THROW(
        check_cholesky_factor_corr("checkCorrCholeskyMatrix", "y", y_vec),
        std::domain_error);
    for (auto& y : y_vec) {
      y << 1, 0, 0, sqrt(0.5), sqrt(0.5), 0, sqrt(0.25), sqrt(0.25), sqrt(0.5);
    }
  }
}
