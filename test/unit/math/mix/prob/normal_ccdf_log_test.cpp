#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/normal_lcdf_tail_test_helpers.hpp>

TEST_F(AgradRev, mathMixScalFun_normal_lccdf) {
  auto f = [](const double mu, const double sigma) {
    return
        [=](const auto& y) { return stan::math::normal_lccdf(y, mu, sigma); };
  };

  stan::test::expect_ad(f(0.0, 1.0), 50.0);
  stan::test::expect_ad(f(0.0, 1.0), 20.0 * stan::math::SQRT_TWO);
  stan::test::expect_ad(f(0.0, 1.0), 5.5);
  stan::test::expect_ad(f(0.0, 1.0), 0.0);
  stan::test::expect_ad(f(0.0, 1.0), -0.15);
  stan::test::expect_ad(f(0.0, 1.0), -1.14);
  stan::test::expect_ad(f(0.0, 1.0), -3.00);
  stan::test::expect_ad(f(0.0, 1.0), -10.00);
  stan::test::expect_ad(f(-1.0, 2.0), -3.50);
  stan::test::expect_ad(f(2.0, 1.0), 3.50);

  // third order autodiff tests can fail at borders of piecewise function
  // stan::test::ad_tolerances tols;
  // tols.grad_hessian_grad_hessian_ = 1e1;
  // stan::test::expect_ad(tols, f(0.0, 1.0), 0.1 * stan::math::SQRT_TWO);
}

namespace {
auto lccdf
    = [](const auto& y) { return stan::math::normal_lccdf(y, 0.0, 1.0); };
constexpr double dir = normal_lcdf_tail_test::orientation::lccdf;
}  // namespace

TEST_F(AgradRev, mathMixScalFun_normal_lccdf_defect_inputs) {
  normal_lcdf_tail_test::expect_ad_at_defect_inputs(lccdf, dir);
}

TEST_F(AgradRev, mathMixScalFun_normal_lccdf_branch_cutoffs) {
  normal_lcdf_tail_test::expect_ad_across_cutoffs(lccdf, dir);
}

TEST_F(AgradRev, mathMixScalFun_normal_lccdf_tail_derivatives) {
  normal_lcdf_tail_test::expect_tail_derivatives(lccdf, dir);
}

TEST_F(AgradRev, mathMixScalFun_normal_lccdf_derivatives_are_finite) {
  normal_lcdf_tail_test::expect_derivatives_finite(lccdf);
}

TEST_F(AgradRev, mathMixScalFun_normal_lccdf_far_tail_gradient) {
  normal_lcdf_tail_test::expect_far_tail_gradient(lccdf, dir);
}

TEST_F(AgradRev, mathMixScalFun_normal_lccdf_branch_accuracy) {
  normal_lcdf_tail_test::expect_branch_accuracy(lccdf, dir);
}
