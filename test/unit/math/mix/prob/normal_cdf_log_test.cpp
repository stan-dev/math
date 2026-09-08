#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/normal_lcdf_tail_test_helpers.hpp>

TEST_F(AgradRev, mathMixScalFun_normal_lcdf) {
  auto f = [](const auto& y, const auto& mu, const auto& sigma) {
    return stan::math::normal_lcdf(y, mu, sigma);
  };

  stan::test::expect_ad(f, -50.0, 0.0, 1.0);
  stan::test::expect_ad(f, -20.0 * stan::math::SQRT_TWO, 0.0, 1.0);
  stan::test::expect_ad(f, -5.5, 0.0, 1.0);
  stan::test::expect_ad(f, 0.0, 0.0, 1.0);
  stan::test::expect_ad(f, 0.15, 0.0, 1.0);
  stan::test::expect_ad(f, 1.14, 0.0, 1.0);
  stan::test::expect_ad(f, 3.00, 0.0, 1.0);
  stan::test::expect_ad(f, 10.00, 0.0, 1.0);
  stan::test::expect_ad(f, 1.50, -1.0, 2.0);
  stan::test::expect_ad(f, 0.50, 2.0, 1.0);
}

namespace normal_lcdf_mix_test {
auto fn = [](const auto& y) { return stan::math::normal_lcdf(y, 0.0, 1.0); };
constexpr double dir = normal_lcdf_tail_test::orientation::lcdf;
}  // namespace normal_lcdf_mix_test

TEST_F(AgradRev, mathMixScalFun_normal_lcdf_defect_inputs) {
  normal_lcdf_tail_test::expect_ad_at_defect_inputs(normal_lcdf_mix_test::fn,
                                                    normal_lcdf_mix_test::dir);
}

TEST_F(AgradRev, mathMixScalFun_normal_lcdf_branch_cutoffs) {
  normal_lcdf_tail_test::expect_ad_across_cutoffs(normal_lcdf_mix_test::fn,
                                                  normal_lcdf_mix_test::dir);
}

TEST_F(AgradRev, mathMixScalFun_normal_lcdf_tail_derivatives) {
  normal_lcdf_tail_test::expect_tail_derivatives(normal_lcdf_mix_test::fn,
                                                 normal_lcdf_mix_test::dir);
}

TEST_F(AgradRev, mathMixScalFun_normal_lcdf_derivatives_are_finite) {
  normal_lcdf_tail_test::expect_derivatives_finite(normal_lcdf_mix_test::fn);
}

TEST_F(AgradRev, mathMixScalFun_normal_lcdf_far_tail_gradient) {
  normal_lcdf_tail_test::expect_far_tail_gradient(normal_lcdf_mix_test::fn,
                                                  normal_lcdf_mix_test::dir);
}

TEST_F(AgradRev, mathMixScalFun_normal_lcdf_branch_accuracy) {
  normal_lcdf_tail_test::expect_branch_accuracy(normal_lcdf_mix_test::fn,
                                                normal_lcdf_mix_test::dir);
}
