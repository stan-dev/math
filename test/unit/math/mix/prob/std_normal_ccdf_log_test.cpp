#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/normal_lcdf_tail_test_helpers.hpp>

TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf) {
  auto f = [](const auto& y) { return stan::math::std_normal_lccdf(y); };

  stan::test::expect_ad(f, 50.0);
  stan::test::expect_ad(f, 20.0 * stan::math::SQRT_TWO);
  stan::test::expect_ad(f, 5.5);
  stan::test::expect_ad(f, 0.0);
  stan::test::expect_ad(f, -0.15);
  stan::test::expect_ad(f, -1.14);
  stan::test::expect_ad(f, -3.00);
  stan::test::expect_ad(f, -10.00);
}

namespace std_normal_lccdf_mix_test {
auto fn = [](const auto& y) { return stan::math::std_normal_lccdf(y); };
constexpr double dir = normal_lcdf_tail_test::orientation::lccdf;
}  // namespace std_normal_lccdf_mix_test

TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf_defect_inputs) {
  normal_lcdf_tail_test::expect_ad_at_defect_inputs(
      std_normal_lccdf_mix_test::fn, std_normal_lccdf_mix_test::dir);
}

TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf_branch_cutoffs) {
  normal_lcdf_tail_test::expect_ad_across_cutoffs(
      std_normal_lccdf_mix_test::fn, std_normal_lccdf_mix_test::dir);
}

TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf_tail_derivatives) {
  normal_lcdf_tail_test::expect_tail_derivatives(
      std_normal_lccdf_mix_test::fn, std_normal_lccdf_mix_test::dir);
}

TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf_derivatives_are_finite) {
  normal_lcdf_tail_test::expect_derivatives_finite(
      std_normal_lccdf_mix_test::fn);
}

TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf_far_tail_gradient) {
  normal_lcdf_tail_test::expect_far_tail_gradient(
      std_normal_lccdf_mix_test::fn, std_normal_lccdf_mix_test::dir);
}

TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf_branch_accuracy) {
  normal_lcdf_tail_test::expect_branch_accuracy(std_normal_lccdf_mix_test::fn,
                                                std_normal_lccdf_mix_test::dir);
}

TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf_varmat) {
  using stan::math::var;
  using stan::math::var_value;
  Eigen::VectorXd y(3);
  y << -1.5, 0.2, 3.0;
  Eigen::Matrix<var, -1, 1> y_mat = y;
  var lp_mat = stan::math::std_normal_lccdf(y_mat);
  lp_mat.grad();
  const Eigen::VectorXd y_adj = y_mat.adj();
  stan::math::set_zero_all_adjoints();

  var_value<Eigen::VectorXd> y_var(y);
  var lp = stan::math::std_normal_lccdf(y_var);
  lp.grad();
  EXPECT_DOUBLE_EQ(lp_mat.val(), lp.val());
  EXPECT_MATRIX_NEAR(y_adj, y_var.adj(), 1e-14);
}
