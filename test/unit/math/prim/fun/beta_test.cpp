#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <test/unit/math/prim/fun/ternary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST(MathFunctions, beta) {
  using stan::math::beta;

  EXPECT_FLOAT_EQ(beta(2.15, 1.71), 0.1936023967178879658641281697269);
  EXPECT_FLOAT_EQ(beta(7.62, 10.15), 0.0000065340071071564116445887286);
}

TEST(MathFunctions, beta_nan) {
  using stan::math::beta;
  using stan::math::INFTY;
  using stan::math::NOT_A_NUMBER;

  EXPECT_TRUE(std::isnan(beta(NOT_A_NUMBER, 2.16)));
  EXPECT_TRUE(std::isnan(beta(1.65, NOT_A_NUMBER)));

  EXPECT_TRUE(std::isnan(beta(INFTY, 2.16)));
  EXPECT_TRUE(std::isnan(beta(1.65, INFTY)));
}

TEST(MathFunctions, beta_vec) {
  auto f
      = [](const auto& x1, const auto& x2) { return stan::math::beta(x1, x2); };

  Eigen::VectorXd in1 = Eigen::VectorXd::Random(6).cwiseAbs();
  Eigen::VectorXd in2 = Eigen::VectorXd::Random(6).cwiseAbs();

  stan::test::binary_scalar_tester(f, in1, in2);
}

TEST(MathFunctions, beta_inc_unnorm) {
  using stan::math::beta;

  EXPECT_FLOAT_EQ(beta(2.15, 1.71, 0.2), 0.0131642393318922);
  EXPECT_FLOAT_EQ(beta(2.15, 1.71, 0.8), 0.161649383778871);
  EXPECT_FLOAT_EQ(beta(7.62, 10.15, 0.9), 6.53400339352810e-6);
}

TEST(MathFunctions, beta_inc_unnorm_vec) {
  auto f = [](const auto& a, const auto& b, const auto& x) {
    return stan::math::beta(a, b, x); };

  Eigen::VectorXd a = Eigen::VectorXd::Random(6).cwiseAbs();
  Eigen::VectorXd b = Eigen::VectorXd::Random(6).cwiseAbs();
  Eigen::VectorXd x = Eigen::VectorXd::Random(6).array().logistic().matrix();

  stan::test::ternary_scalar_tester(f, a, b, x);
}
