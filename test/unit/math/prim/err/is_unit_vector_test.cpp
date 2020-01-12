#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>

TEST(ErrorHandlingMatrix, isUnitVector) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y(2);

  y << sqrt(0.5), sqrt(0.5);
  EXPECT_TRUE(stan::math::is_unit_vector(y));

  y[1] = 0;
  EXPECT_FALSE(stan::math::is_unit_vector(y));
}

TEST(ErrorHandlingMatrix, isUnitVector_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y(2);
  double nan = std::numeric_limits<double>::quiet_NaN();

  y << nan, sqrt(0.5);
  EXPECT_FALSE(stan::math::is_unit_vector(y));

  y << sqrt(0.5), nan;
  EXPECT_FALSE(stan::math::is_unit_vector(y));

  y << nan, nan;
  EXPECT_FALSE(stan::math::is_unit_vector(y));
}

TEST(ErrorHandlingMatrix, isUnitVector_0_size) {
  using stan::math::is_unit_vector;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y(0, 1);

  EXPECT_FALSE(is_unit_vector(y));
}
