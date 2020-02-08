#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <limits>

TEST(ErrorHandlingMatrix, checkMultiplicablePositiveMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;

  y.resize(3, 3);
  x.resize(3, 3);
  EXPECT_NO_THROW(stan::math::check_multiplicable_positive(
      "checkMultiplicablePositive", "x", x, "y", y));
  x.resize(3, 2);
  y.resize(2, 4);
  EXPECT_NO_THROW(stan::math::check_multiplicable_positive(
      "checkMultiplicablePositive", "x", x, "y", y));

  y.resize(1, 2);
  EXPECT_THROW(stan::math::check_multiplicable_positive(
                   "checkMultiplicablePositive", "x", x, "y", y),
               std::invalid_argument);

  x.resize(2, 2);
  EXPECT_THROW(stan::math::check_multiplicable_positive(
                   "checkMultiplicablePositive", "x", x, "y", y),
               std::invalid_argument);
}

TEST(ErrorHandlingMatrix, checkMultiplicablePositiveMatrix_0) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;

  x.resize(3, 0);
  y.resize(0, 3);
  EXPECT_THROW(stan::math::check_multiplicable_positive(
                   "checkMultiplicablePositive", "x", x, "y", y),
               std::invalid_argument);

  x.resize(0, 4);
  y.resize(4, 3);
  EXPECT_THROW(stan::math::check_multiplicable_positive(
                   "checkMultiplicablePositive", "x", x, "y", y),
               std::invalid_argument);

  x.resize(3, 4);
  y.resize(4, 0);
  EXPECT_THROW(stan::math::check_multiplicable_positive(
                   "checkMultiplicablePositive", "x", x, "y", y),
               std::invalid_argument);
}

TEST(ErrorHandlingMatrix, checkMultiplicablePositiveMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(3, 3);
  x.resize(3, 3);
  y << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  x << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  EXPECT_NO_THROW(stan::math::check_multiplicable_positive(
      "checkMultiplicablePositive", "x", x, "y", y));
  x.resize(3, 2);
  y.resize(2, 4);
  y << nan, nan, nan, nan, nan, nan, nan, nan;
  x << nan, nan, nan, nan, nan, nan;
  EXPECT_NO_THROW(stan::math::check_multiplicable_positive(
      "checkMultiplicablePositive", "x", x, "y", y));

  y.resize(1, 2);
  y << nan, nan;
  EXPECT_THROW(stan::math::check_multiplicable_positive(
                   "checkMultiplicablePositive", "x", x, "y", y),
               std::invalid_argument);

  x.resize(2, 2);
  x << nan, nan, nan, nan;
  EXPECT_THROW(stan::math::check_multiplicable_positive(
                   "checkMultiplicablePositive", "x", x, "y", y),
               std::invalid_argument);
}
