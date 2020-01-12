#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, checkCovMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(3, 3);
  y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  EXPECT_NO_THROW(stan::math::check_cov_matrix("checkCovMatrix", "y", y));

  y << 1, 2, 3, 2, 1, 2, 3, 2, 1;
  EXPECT_THROW(stan::math::check_cov_matrix("checkCovMatrix", "y", y),
               std::domain_error);
}

TEST(ErrorHandlingMatrix, checkCovMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(3, 3);
  y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  EXPECT_NO_THROW(stan::math::check_cov_matrix("checkCovMatrix", "y", y));

  for (int i = 0; i < y.size(); i++) {
    y.resize(3, 3);
    y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
    y(i) = nan;
    EXPECT_THROW(stan::math::check_cov_matrix("checkCovMatrix", "y", y),
                 std::domain_error);
    y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  }
}
