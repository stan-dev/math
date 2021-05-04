#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, isLDLTFactor_nan) {
  using stan::math::is_ldlt_factor;

  double nan = std::numeric_limits<double>::quiet_NaN();

  {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x(2, 2);
    x << 2, 1, 1, 2;
    auto ldlt_x = stan::math::make_ldlt_factor<Eigen::MatrixXd>(x);
    EXPECT_TRUE(is_ldlt_factor(ldlt_x));
  }

  {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x(2, 2);
    x << nan, 1, 1, 3;
    auto ldlt_x = stan::math::make_ldlt_factor<Eigen::MatrixXd>(x);
    EXPECT_FALSE(is_ldlt_factor(ldlt_x));
  }

  {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x(2, 2);
    x << 3, nan, 1, 3;
    auto ldlt_x = stan::math::make_ldlt_factor<Eigen::MatrixXd>(x);
    EXPECT_TRUE(is_ldlt_factor(ldlt_x));
  }

  {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x(2, 2);
    x << 3, 1, nan, 3;
    auto ldlt_x = stan::math::make_ldlt_factor<Eigen::MatrixXd>(x);
    EXPECT_FALSE(is_ldlt_factor(ldlt_x));
  }

  {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x(2, 2);
    x << 3, 1, 1, nan;
    auto ldlt_x = stan::math::make_ldlt_factor<Eigen::MatrixXd>(x);
    EXPECT_FALSE(is_ldlt_factor(ldlt_x));
  }
}
