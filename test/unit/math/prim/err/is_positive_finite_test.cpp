#include <gtest/gtest.h>
#include <stan/math/prim/err/is_positive_finite.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <vector>

using stan::math::is_positive_finite;

TEST(ErrorHandlingScalar, isPositiveFinite) {
  EXPECT_FALSE(is_positive_finite(stan::math::NEGATIVE_INFTY));
  EXPECT_FALSE(is_positive_finite(-3.0));
  EXPECT_FALSE(is_positive_finite(-3));
  EXPECT_FALSE(is_positive_finite(-0.0));
  EXPECT_FALSE(is_positive_finite(0.0));
  EXPECT_FALSE(is_positive_finite(0u));
  EXPECT_FALSE(is_positive_finite((size_t)0));
  EXPECT_FALSE(is_positive_finite(0));
  EXPECT_TRUE(is_positive_finite(3.0));
  EXPECT_TRUE(is_positive_finite(3));
  EXPECT_FALSE(is_positive_finite(stan::math::INFTY));
  EXPECT_FALSE(is_positive_finite(stan::math::NOT_A_NUMBER));
}

TEST(ErrorHandlingScalar, isPositiveFiniteVectorization) {
  Eigen::MatrixXd m = Eigen::MatrixXd::Constant(3, 2, 1);
  EXPECT_TRUE(is_positive_finite(std::vector<Eigen::MatrixXd>{m, m, m}));
  EXPECT_TRUE(is_positive_finite(m));
  Eigen::MatrixXd m2 = m;
  m2(1, 1) = -0.5;
  EXPECT_FALSE(is_positive_finite(std::vector<Eigen::MatrixXd>{m, m2, m}));
  EXPECT_FALSE(is_positive_finite(m2));
  m2(1, 1) = stan::math::NOT_A_NUMBER;
  EXPECT_FALSE(is_positive_finite(std::vector<Eigen::MatrixXd>{m, m2, m}));
  EXPECT_FALSE(is_positive_finite(m2));
}
