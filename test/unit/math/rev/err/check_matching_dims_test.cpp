#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>

TEST(ErrorHandlingMatrix, checkMatchingDims_double_var) {
  using stan::math::check_matching_dims;

  std::vector<Eigen::Matrix<double, -1, 1>> x;
  std::vector<Eigen::Matrix<stan::math::var, -1, 1>> y;
  x = std::vector<Eigen::Matrix<double, -1, 1>>(
      5, Eigen::Matrix<double, -1, 1>(4));
  y = std::vector<Eigen::Matrix<stan::math::var, -1, 1>>(
      5, Eigen::Matrix<stan::math::var, -1, 1>(4));

  EXPECT_NO_THROW(check_matching_dims("checkMatchingDims", "x", x, "y", y));

  x = std::vector<Eigen::Matrix<double, -1, 1>>(
      5, Eigen::Matrix<double, -1, 1>(6));
  y = std::vector<Eigen::Matrix<stan::math::var, -1, 1>>(
      5, Eigen::Matrix<stan::math::var, -1, 1>(4));

  EXPECT_THROW(check_matching_dims("checkMatchingDims", "x", x, "y", y),
               std::invalid_argument);
}

TEST(ErrorHandlingMatrix, checkMatchingDims_var) {
  using stan::math::check_matching_dims;

  std::vector<Eigen::Matrix<stan::math::var, -1, 1>> x;
  std::vector<Eigen::Matrix<stan::math::var, -1, 1>> y;
  x = std::vector<Eigen::Matrix<stan::math::var, -1, 1>>(
      5, Eigen::Matrix<stan::math::var, -1, 1>(4));
  y = std::vector<Eigen::Matrix<stan::math::var, -1, 1>>(
      5, Eigen::Matrix<stan::math::var, -1, 1>(4));

  EXPECT_NO_THROW(check_matching_dims("checkMatchingDims", "x", x, "y", y));

  x = std::vector<Eigen::Matrix<stan::math::var, -1, 1>>(
      5, Eigen::Matrix<stan::math::var, -1, 1>(6));
  y = std::vector<Eigen::Matrix<stan::math::var, -1, 1>>(
      5, Eigen::Matrix<stan::math::var, -1, 1>(4));

  EXPECT_THROW(check_matching_dims("checkMatchingDims", "x", x, "y", y),
               std::invalid_argument);
}
