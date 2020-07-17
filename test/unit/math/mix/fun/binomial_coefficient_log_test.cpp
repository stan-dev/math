#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, binomialCoefficientLog) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::binomial_coefficient_log(x1, x2);
  };
  stan::test::expect_ad(f, 3, 2);
  stan::test::expect_ad(f, 24.0, 12.0);
  stan::test::expect_ad(f, 1.0, 0.0);
  stan::test::expect_ad(f, 0.0, 1.0);
  stan::test::expect_ad(f, -0.3, 0.5);

  stan::test::expect_common_nonzero_binary(f);
}

TEST(mathMixScalFun, binomialCoefficientLog_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::binomial_coefficient_log;
    return binomial_coefficient_log(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 3.0, 1.5;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}
