#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>

TEST(ErrorHandlingMatrix, checkUnitVector) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(2);
  y_vec << sqrt(0.5), sqrt(0.5);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};

  EXPECT_NO_THROW(stan::math::check_unit_vector("checkUnitVector", "y", y));
  for (auto& y_i : y) {
    y_i[1] = 0;
  }
  EXPECT_THROW(stan::math::check_unit_vector("checkUnitVector", "y", y),
               std::domain_error);
}

TEST(ErrorHandlingMatrix, checkUnitVector_nan) {
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(2);
  y_vec << nan, sqrt(0.5);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};

  EXPECT_THROW(stan::math::check_unit_vector("checkUnitVector", "y", y),
               std::domain_error);
  for (auto& y_i : y) {
    y_i << sqrt(0.5), nan;
  }
  EXPECT_THROW(stan::math::check_unit_vector("checkUnitVector", "y", y),
               std::domain_error);
  for (auto& y_i : y) {
    y_i << nan, nan;
  }
  EXPECT_THROW(stan::math::check_unit_vector("checkUnitVector", "y", y),
               std::domain_error);
}

TEST(ErrorHandlingMatrix, checkUnitVector_0_size) {
  using stan::math::check_unit_vector;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(0, 1);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};

  EXPECT_THROW_MSG(check_unit_vector("checkUnitVector", "y", y),
                   std::invalid_argument,
                   "y[1] has size 0, but must have a non-zero size");
}
