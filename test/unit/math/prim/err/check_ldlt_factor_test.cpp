#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, CheckLDLTFactor_nan) {
  using stan::math::check_ldlt_factor;

  double nan = std::numeric_limits<double>::quiet_NaN();
  {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x(2, 2);
    x << 2, 1, 1, 2;
    auto ldlt_x = stan::math::make_ldlt_factor<Eigen::MatrixXd>(x);
    EXPECT_NO_THROW(
        check_ldlt_factor("checkLDLTFactorMatrix", "ldlt_x", ldlt_x));
  }

  {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x(2, 2);
    x << nan, 1, 1, 3;
    auto ldlt_x = stan::math::make_ldlt_factor<Eigen::MatrixXd>(x);
    EXPECT_THROW(check_ldlt_factor("checkLDLTFactorMatrix", "ldlt_x", ldlt_x),
                 std::domain_error);
  }

  {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x(2, 2);
    x << 3, nan, 1, 3;
    auto ldlt_x = stan::math::make_ldlt_factor<Eigen::MatrixXd>(x);
    EXPECT_NO_THROW(
        check_ldlt_factor("checkLDLTFactorMatrix", "ldlt_x", ldlt_x));
  }

  {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x(2, 2);
    x << 3, 1, nan, 3;
    auto ldlt_x = stan::math::make_ldlt_factor<Eigen::MatrixXd>(x);
    EXPECT_THROW(check_ldlt_factor("checkLDLTFactorMatrix", "ldlt_x", ldlt_x),
                 std::domain_error);
  }

  {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x(2, 2);
    x << 3, 1, 1, nan;
    auto ldlt_x = stan::math::make_ldlt_factor<Eigen::MatrixXd>(x);
    EXPECT_THROW(check_ldlt_factor("checkLDLTFactorMatrix", "ldlt_x", ldlt_x),
                 std::domain_error);
  }
}
