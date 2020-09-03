#include <gtest/gtest.h>
#include <stan/math/prim/err/is_not_nan.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <vector>

TEST(ErrorHandlingScalar, isNotNan) {
  using stan::math::is_not_nan;

  EXPECT_TRUE(is_not_nan(stan::math::NEGATIVE_INFTY));
  EXPECT_TRUE(is_not_nan(-3.0));
  EXPECT_TRUE(is_not_nan(-3));
  EXPECT_TRUE(is_not_nan(-0.0));
  EXPECT_TRUE(is_not_nan(0.0));
  EXPECT_TRUE(is_not_nan(0u));
  EXPECT_TRUE(is_not_nan((size_t)0));
  EXPECT_TRUE(is_not_nan(0));
  EXPECT_TRUE(is_not_nan(3.0));
  EXPECT_TRUE(is_not_nan(3));
  EXPECT_TRUE(is_not_nan(stan::math::INFTY));
  EXPECT_TRUE(is_not_nan(-stan::math::INFTY));
  EXPECT_FALSE(is_not_nan(stan::math::NOT_A_NUMBER));
}

TEST(ErrorHandlingScalar, isNotNanEigen) {
  using stan::math::is_not_nan;

  Eigen::MatrixXd m = Eigen::MatrixXd::Constant(3, 2, 1);
  EXPECT_TRUE(is_not_nan(m));
  Eigen::MatrixXd m2 = m;
  m2(1, 1) = stan::math::NOT_A_NUMBER;
  EXPECT_FALSE(is_not_nan(m2));
}

TEST(ErrorHandlingScalar, isNotNanVectorization) {
  using stan::math::is_not_nan;

  Eigen::MatrixXd m = Eigen::MatrixXd::Constant(3, 2, 1);
  EXPECT_TRUE(is_not_nan(std::vector<Eigen::MatrixXd>{m, m, m}));
  Eigen::MatrixXd m2 = m;
  m2(1, 1) = stan::math::NOT_A_NUMBER;
  EXPECT_FALSE(is_not_nan(std::vector<Eigen::MatrixXd>{m, m2, m}));
}
