#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, hypergeometric_pFq_values) {
  using Eigen::VectorXd;
  using stan::math::hypergeometric_pFq;

  VectorXd a(2);
  VectorXd b(2);
  a << 4, 4;
  b << 5, 5;
  double z = 2;

  EXPECT_FLOAT_EQ(3.8420514314107791, hypergeometric_pFq(a, b, z));

  a << 6, 4;
  b << 3, 1;
  z = 2;

  EXPECT_FLOAT_EQ(848.09943891059574, hypergeometric_pFq(a, b, z));
}

TEST(MathFunctions, hypergeometric_pFq_errors) {
  using Eigen::VectorXd;
  using stan::math::hypergeometric_pFq;
  using stan::math::INFTY;
  using stan::math::NOT_A_NUMBER;

  VectorXd a(2);
  VectorXd b(2);
  a << 4, INFTY;
  b << 5, 5;
  double z = 2;

  EXPECT_THROW(hypergeometric_pFq(a, b, z), std::domain_error);

  a << 6, 4;
  b << NOT_A_NUMBER, 1;

  EXPECT_THROW(hypergeometric_pFq(a, b, z), std::domain_error);
}
