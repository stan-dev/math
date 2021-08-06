#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, checkCovMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_mat;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y_vec{
      y_mat, y_mat, y_mat};

  for (auto& y : y_vec) {
    y.resize(3, 3);
    y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  }
  EXPECT_NO_THROW(stan::math::check_cov_matrix("checkCovMatrix", "y", y_vec));

  for (auto& y : y_vec) {
    y << 1, 2, 3, 2, 1, 2, 3, 2, 1;
  }
  EXPECT_THROW(stan::math::check_cov_matrix("checkCovMatrix", "y", y_vec),
               std::domain_error);
}

TEST(ErrorHandlingMatrix, checkCovMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_mat;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y_vec{
      y_mat, y_mat, y_mat};
  double nan = std::numeric_limits<double>::quiet_NaN();

  for (auto& y : y_vec) {
    y.resize(3, 3);
    y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  }
  EXPECT_NO_THROW(stan::math::check_cov_matrix("checkCovMatrix", "y", y_vec));

  for (int i = 0; i < y_vec[0].size(); i++) {
    for (auto& y : y_vec) {
      y.resize(3, 3);
      y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
      y(i) = nan;
    }
    EXPECT_THROW(stan::math::check_cov_matrix("checkCovMatrix", "y", y_vec),
                 std::domain_error);
    for (auto& y : y_vec) {
      y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
    }
  }
}
