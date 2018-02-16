#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>

TEST(ErrorHandlingMatrix, checkUnitVector) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y(2);
  y << sqrt(0.5), sqrt(0.5);

  EXPECT_NO_THROW(stan::math::check_unit_vector("checkUnitVector", "y", y));

  y[1] = 0;
  EXPECT_THROW(stan::math::check_unit_vector("checkUnitVector", "y", y),
               std::domain_error);
}

TEST(ErrorHandlingMatrix, checkUnitVector_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y(2);
  double nan = std::numeric_limits<double>::quiet_NaN();

  y << nan, sqrt(0.5);
  EXPECT_THROW(stan::math::check_unit_vector("checkUnitVector", "y", y),
               std::domain_error);
  y << sqrt(0.5), nan;
  EXPECT_THROW(stan::math::check_unit_vector("checkUnitVector", "y", y),
               std::domain_error);
  y << nan, nan;
  EXPECT_THROW(stan::math::check_unit_vector("checkUnitVector", "y", y),
               std::domain_error);
}

TEST(ErrorHandlingMatrix, checkUnitVector_0_size) {
  using stan::math::check_unit_vector;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y(0, 1);

  EXPECT_THROW_MSG(check_unit_vector("checkUnitVector", "y", y),
                   std::invalid_argument,
                   "y has size 0, but must have a non-zero size");
}
