#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, logDiffExp) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::log_diff_exp(x1, x2);
  };
  // stan::test::expect_common_nonzero_binary(f);
  stan::test::expect_ad(f, 1.2, 0.5);
  stan::test::expect_ad(f, 5.0, 2.0);
  stan::test::expect_ad(f, 9.0, 6.0);
  stan::test::expect_ad(f, 0.1, 0.0);
  stan::test::expect_ad(f, 3.01, 2.0);
  stan::test::expect_ad(f, 2.0, 1.0);

  // hessians near zero have high relative error but small absolute error
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 3.0;
  tols.hessian_fvar_hessian_ = 3.0;
  tols.grad_hessian_hessian_ = 3.0;
  tols.grad_hessian_grad_hessian_ = 3.0;
  stan::test::expect_ad(tols, f, -2.0, -15.0);
  stan::test::expect_ad(tols, f, 2.0, -15.0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);
}

TEST(mathMixScalFun, logDiffExp_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::log_diff_exp;
    return log_diff_exp(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}
