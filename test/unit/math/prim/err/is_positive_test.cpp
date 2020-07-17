#include <gtest/gtest.h>
#include <stan/math/prim/err/is_positive.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <vector>

using stan::math::is_positive;

TEST(ErrorHandlingScalar, isPositive) {
  EXPECT_FALSE(is_positive(stan::math::NEGATIVE_INFTY));
  EXPECT_FALSE(is_positive(-3.0));
  EXPECT_FALSE(is_positive(-3));
  EXPECT_FALSE(is_positive(-0.0));
  EXPECT_FALSE(is_positive(0.0));
  EXPECT_FALSE(is_positive(0u));
  EXPECT_FALSE(is_positive((size_t)0));
  EXPECT_FALSE(is_positive(0));
  EXPECT_TRUE(is_positive(3.0));
  EXPECT_TRUE(is_positive(3));
  EXPECT_TRUE(is_positive(stan::math::INFTY));
  EXPECT_FALSE(is_positive(stan::math::NOT_A_NUMBER));
}

TEST(ErrorHandlingScalar, isPositiveVectorization) {
  Eigen::MatrixXd m = Eigen::MatrixXd::Constant(3, 2, 1);
  EXPECT_TRUE(is_positive(std::vector<Eigen::MatrixXd>{m, m, m}));
  Eigen::MatrixXd m2 = m;
  m2(1, 1) = -0.5;
  EXPECT_FALSE(is_positive(std::vector<Eigen::MatrixXd>{m, m2, m}));
}
