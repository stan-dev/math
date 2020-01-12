#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, checkCovCholeskyMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  using stan::math::check_cholesky_factor;

  y.resize(1, 1);
  y << 1;
  EXPECT_NO_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y));

  y.resize(3, 3);
  y << 1, 0, 0, 1, 1, 0, 1, 1, 1;
  EXPECT_NO_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y));

  // not positive
  y.resize(1, 1);
  y << -1;
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y),
               std::domain_error);

  // not lower triangular
  y.resize(3, 3);
  y << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y),
               std::domain_error);

  // not positive
  y.resize(3, 3);
  y << 1, 0, 0, 2, -1, 0, 1, 2, 3;
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y),
               std::domain_error);

  // not rectangular
  y.resize(2, 3);
  y << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y),
               std::domain_error);
  y.resize(3, 2);
  y << 1, 0, 2, 3, 4, 5;
  EXPECT_NO_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y));
  y(0, 1) = 1.5;
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y),
               std::domain_error);
}

TEST(ErrorHandlingMatrix, checkCovCholeskyMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();

  using stan::math::check_cholesky_factor;

  y.resize(1, 1);
  y << nan;
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y),
               std::domain_error);

  y.resize(3, 3);
  y << nan, 0, 0, nan, nan, 0, nan, nan, nan;
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y),
               std::domain_error);

  // not positive
  y.resize(1, 1);
  y << nan;
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y),
               std::domain_error);

  // not lower triangular
  y.resize(3, 3);
  y << 1, nan, 3, 4, 5, nan, 7, 8, 9;
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y),
               std::domain_error);

  // not positive
  y.resize(3, 3);
  y << 1, 0, 0, 2, nan, 0, 1, 2, 3;
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y),
               std::domain_error);

  // not rectangular
  y.resize(2, 3);
  y << 1, 2, nan, nan, 5, 6;
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y),
               std::domain_error);
  y.resize(3, 2);
  y << 1, 0, 2, nan, 4, 5;
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y),
               std::domain_error);
  y(0, 1) = nan;
  EXPECT_THROW(check_cholesky_factor("checkCovCholeskyMatrix", "y", y),
               std::domain_error);
}
