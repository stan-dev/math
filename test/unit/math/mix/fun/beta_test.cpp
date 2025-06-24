#include <test/unit/math/test_ad.hpp>
/*
TEST(mathMixScalFun, beta) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::beta;
    return beta(x1, x2);
  };

  stan::test::expect_common_nonzero_binary(f);

  stan::test::expect_ad(f, 0.5, 3.3);
  stan::test::expect_ad(f, 3.4, 0.9);
  stan::test::expect_ad(f, 5.2, 6.7);
  stan::test::expect_ad(f, 7.5, 1.8);

  Eigen::VectorXd in1(3);
  in1 << 0.5, 3.4, 5.2;
  Eigen::VectorXd in2(3);
  in2 << 3.3, 0.9, 6.7;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}

TEST(mathMixScalFun, beta_varmat) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::beta;
    return beta(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 0.5, 3.4, 5.2;
  Eigen::VectorXd in2(3);
  in2 << 3.3, 0.9, 6.7;
  stan::test::expect_ad_matvar(f, in1, in2);
  stan::test::expect_ad_matvar(f, in1(0), in2);
  stan::test::expect_ad_matvar(f, in1, in2(0));
}
*/
TEST(mathMixScalFun, beta_inc_unnorm) {
  auto f = [](const auto& a, const auto& b, const auto& x) {
    using stan::math::beta;
    using stan::math::exp;
    using stan::math::inv_logit;

    return beta(exp(a), exp(b), inv_logit(x));
  };

  stan::test::expect_ad(f, 0.5, 1.2, 0.7);
  stan::test::expect_ad(f, 1.1, 3.4, 0.2);
  stan::test::expect_ad(f, 2.3, 0.9, 0.4);
}
