#include <gtest/gtest.h>
#include <stan/math/prim/err/is_scal_finite.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <vector>

TEST(ErrorHandlingScalar, isScalFinite) {
  using stan::math::is_scal_finite;
  EXPECT_FALSE(is_scal_finite(stan::math::NEGATIVE_INFTY));
  EXPECT_TRUE(is_scal_finite(-3.0));
  EXPECT_TRUE(is_scal_finite(-3));
  EXPECT_TRUE(is_scal_finite(-0.0));
  EXPECT_TRUE(is_scal_finite(0.0));
  EXPECT_TRUE(is_scal_finite(0u));
  EXPECT_TRUE(is_scal_finite((size_t)0));
  EXPECT_TRUE(is_scal_finite(0));
  EXPECT_TRUE(is_scal_finite(3.0));
  EXPECT_TRUE(is_scal_finite(3));
  EXPECT_FALSE(is_scal_finite(stan::math::INFTY));
  EXPECT_FALSE(is_scal_finite(-stan::math::INFTY));
  EXPECT_FALSE(is_scal_finite(stan::math::NOT_A_NUMBER));
}

TEST(ErrorHandlingScalar, isScalFiniteEigen) {
  using stan::math::is_scal_finite;
  Eigen::MatrixXd m = Eigen::MatrixXd::Constant(3, 2, 1);
  EXPECT_TRUE(is_scal_finite(m));
  Eigen::MatrixXd m2 = m;
  m2(1, 1) = stan::math::NOT_A_NUMBER;
  EXPECT_FALSE(is_scal_finite(m2));
}

TEST(ErrorHandlingScalar, isScalFiniteVectorization) {
  using stan::math::is_scal_finite;
  Eigen::MatrixXd m = Eigen::MatrixXd::Constant(3, 2, 1);
  EXPECT_TRUE(is_scal_finite(std::vector<Eigen::MatrixXd>{m, m, m}));
  Eigen::MatrixXd m2 = m;
  m2(1, 1) = stan::math::NOT_A_NUMBER;
  EXPECT_FALSE(is_scal_finite(std::vector<Eigen::MatrixXd>{m, m2, m}));
}
