#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, checkCovCholeskyMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_mat;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y_vec{
      y_mat, y_mat, y_mat};

  using stan::math::check_cholesky_factor;
  for (auto& y : y_vec) {
    y.resize(1, 1);
    y << 1;
  }
  EXPECT_NO_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec));

  for (auto& y : y_vec) {
    y.resize(3, 3);
    y << 1, 0, 0, 1, 1, 0, 1, 1, 1;
  }
  EXPECT_NO_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec));

  // not positive
  for (auto& y : y_vec) {
    y.resize(1, 1);
    y << -1;
  }
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec),
               std::domain_error);

  // not lower triangular
  for (auto& y : y_vec) {
    y.resize(3, 3);
    y << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  }
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec),
               std::domain_error);

  // not positive
  for (auto& y : y_vec) {
    y.resize(3, 3);
    y << 1, 0, 0, 2, -1, 0, 1, 2, 3;
  }
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec),
               std::domain_error);

  // not rectangular
  for (auto& y : y_vec) {
    y.resize(2, 3);
    y << 1, 2, 3, 4, 5, 6;
  }
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec),
               std::domain_error);
  for (auto& y : y_vec) {
    y.resize(3, 2);
    y << 1, 0, 2, 3, 4, 5;
  }
  EXPECT_NO_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec));
  for (auto& y : y_vec) {
    y(0, 1) = 1.5;
  }
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec),
               std::domain_error);
}

TEST(ErrorHandlingMatrix, checkCovCholeskyMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_mat;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y_vec{
      y_mat, y_mat, y_mat};
  double nan = std::numeric_limits<double>::quiet_NaN();

  using stan::math::check_cholesky_factor;
  for (auto& y : y_vec) {
    y.resize(1, 1);
    y << nan;
  }
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec),
               std::domain_error);

  for (auto& y : y_vec) {
    y.resize(3, 3);
    y << nan, 0, 0, nan, nan, 0, nan, nan, nan;
  }
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec),
               std::domain_error);

  // not positive
  for (auto& y : y_vec) {
    y.resize(1, 1);
    y << nan;
  }
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec),
               std::domain_error);

  // not lower triangular
  for (auto& y : y_vec) {
    y.resize(3, 3);
    y << 1, nan, 3, 4, 5, nan, 7, 8, 9;
  }
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec),
               std::domain_error);

  // not positive
  for (auto& y : y_vec) {
    y.resize(3, 3);
    y << 1, 0, 0, 2, nan, 0, 1, 2, 3;
  }
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec),
               std::domain_error);

  // not rectangular
  for (auto& y : y_vec) {
    y.resize(2, 3);
    y << 1, 2, nan, nan, 5, 6;
  }
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec),
               std::domain_error);
  for (auto& y : y_vec) {
    y.resize(3, 2);
    y << 1, 0, 2, nan, 4, 5;
  }
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec),
               std::domain_error);
  for (auto& y : y_vec) {
    y(0, 1) = nan;
  }
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y_vec),
               std::domain_error);
}
