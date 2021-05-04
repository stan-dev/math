#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, isCorrCholeskyMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  using stan::math::is_cholesky_factor_corr;
  using std::sqrt;

  y.resize(1, 1);
  y << 1;
  EXPECT_TRUE(is_cholesky_factor_corr(y));

  y.resize(3, 3);
  y << 1, 0, 0, sqrt(0.5), sqrt(0.5), 0, sqrt(0.25), sqrt(0.25), sqrt(0.5);
  EXPECT_TRUE(is_cholesky_factor_corr(y));

  // not positive
  y.resize(1, 1);
  y << -1;
  EXPECT_FALSE(is_cholesky_factor_corr(y));

  // not lower triangular
  y.resize(3, 3);
  y << 1, 2, 3, 0, 5, 6, 0, 0, 9;
  EXPECT_FALSE(is_cholesky_factor_corr(y));

  // not positive
  y.resize(3, 3);
  y << 1, 0, 0, 2, -1, 0, 1, 2, 3;
  EXPECT_FALSE(is_cholesky_factor_corr(y));

  // not rectangular
  y.resize(2, 3);
  y << 1, 2, 3, 4, 5, 6;
  EXPECT_FALSE(is_cholesky_factor_corr(y));
  y.resize(3, 2);
  y << 1, 0, 2, 3, 4, 5;
  EXPECT_FALSE(is_cholesky_factor_corr(y));

  // not unit vectors
  y.resize(3, 3);
  y << 1, 0, 0, 1, 1, 0, 1, 1, 1;
  EXPECT_FALSE(is_cholesky_factor_corr(y));
}

TEST(ErrorHandlingMatrix, isCorrCholeskyMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();

  using stan::math::is_cholesky_factor_corr;
  using std::sqrt;

  y.resize(1, 1);
  y << nan;
  EXPECT_FALSE(is_cholesky_factor_corr(y));

  y.resize(3, 3);
  y << 1, 0, 0, sqrt(0.5), sqrt(0.5), 0, sqrt(0.25), sqrt(0.25), sqrt(0.5);
  EXPECT_TRUE(is_cholesky_factor_corr(y));

  for (int i = 0; i < y.size(); i++) {
    y(i) = nan;
    EXPECT_FALSE(is_cholesky_factor_corr(y));
    y << 1, 0, 0, sqrt(0.5), sqrt(0.5), 0, sqrt(0.25), sqrt(0.25), sqrt(0.5);
  }
}
