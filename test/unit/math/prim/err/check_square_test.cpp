#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, checkSquareMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(3, 3);
  EXPECT_NO_THROW(stan::math::check_square("checkSquareMatrix", "y", y));

  y.resize(3, 2);
  EXPECT_THROW(stan::math::check_square("checkSquareMatrix", "y", y),
               std::invalid_argument);
}

TEST(ErrorHandlingMatrix, checkSquareMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(3, 3);
  y << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  EXPECT_NO_THROW(stan::math::check_square("checkSquareMatrix", "y", y));

  y.resize(3, 2);
  y << nan, nan, nan, nan, nan, nan;
  EXPECT_THROW(stan::math::check_square("checkSquareMatrix", "y", y),
               std::invalid_argument);
}

TEST(ErrorHandlingMatrix, checkSquareMatrix_0x0) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(0, 0);
  EXPECT_NO_THROW(stan::math::check_square("checkSquareMatrix", "y", y));
}

TEST(ErrorHandlingMatrix, checkSquareMatrix_0_size) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(0, 10);
  EXPECT_THROW(stan::math::check_square("checkSquareMatrix", "y", y),
               std::invalid_argument);

  y.resize(10, 0);
  EXPECT_THROW(stan::math::check_square("checkSquareMatrix", "y", y),
               std::invalid_argument);
}
