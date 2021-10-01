#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(MathMixMatFun, logSumExp_signed) {
  auto f = [](const auto x2) {
    return
        [=](const auto& x1) { return stan::math::log_sum_exp_signed(x1, x2); };
  };

  Eigen::VectorXd x0(0);
  Eigen::VectorXi x0_signs(0);
  stan::test::expect_ad(f(x0_signs), x0);
  stan::test::expect_ad_matvar(f(x0_signs), x0);

  Eigen::VectorXd x1(1);
  x1 << 0;
  Eigen::VectorXi x1_signs(1);
  x1_signs << 1;
  stan::test::expect_ad(f(x1_signs), x1);
  stan::test::expect_ad_matvar(f(x1_signs), x1);

  Eigen::VectorXd x2(2);
  x2 << 5, 2;
  Eigen::VectorXi x2_signs(2);
  x2_signs << 1, -1;
  stan::test::expect_ad(f(x2_signs), x2);
  stan::test::expect_ad_matvar(f(x2_signs), x2);

  Eigen::VectorXd x2b(2);
  x2b << 4.9, -std::numeric_limits<double>::infinity();
  Eigen::VectorXi x2b_signs(2);
  x2b_signs << -1, 1;
  stan::test::expect_ad(f(x2b_signs), x2b);
  stan::test::expect_ad_matvar(f(x2b_signs), x2b);

  Eigen::VectorXd x2c(2);
  x2c << std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity();
  Eigen::VectorXi x2c_signs(2);
  x2c_signs << 1, -1;
  stan::test::expect_ad(f(x2c_signs), x2c);
  stan::test::expect_ad_matvar(f(x2c_signs), x2c);

  Eigen::VectorXd x4(4);
  x4 << 1, 2, 3, 4;
  Eigen::VectorXi x4_signs(4);
  x4 << -1, 1, -1, 1;
  stan::test::expect_ad(f(x4_signs), x4);
  stan::test::expect_ad_matvar(f(x4_signs), x4);

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  std::vector<Eigen::VectorXd> inputs{x1, x2, x2b, x2c, x4};
  std::vector<Eigen::VectorXi> signs{x1_signs, x2_signs, x2b_signs, x2c_signs,
                                     x4_signs};

  for (int i = 0; i < 5; i++) {
    stan::test::expect_ad(tols, f(signs[i]), inputs[i]);
    stan::test::expect_ad_matvar(tols, f(signs[i]), inputs[i]);
    Eigen::RowVectorXd rx = inputs[i];
    Eigen::RowVectorXi rsigns = signs[i];
    stan::test::expect_ad(tols, f(rsigns), rx);
    stan::test::expect_ad_matvar(tols, f(rsigns), rx);
    std::vector<double> stx = std::vector<double>(
        inputs[i].data(), inputs[i].data() + inputs[i].size());
    std::vector<int> stsigns = std::vector<int>(
        inputs[i].data(), inputs[i].data() + inputs[i].size());
    stan::test::expect_ad(tols, f(stsigns), stx);
  }

  Eigen::MatrixXd x23(2, 2);
  x23 << 1, 2, 3, 4;
  Eigen::MatrixXi x23_signs(2, 2);
  x23_signs << 1, -1, 1, -1;
  stan::test::expect_ad(f(x23_signs), x23);
  stan::test::expect_ad_matvar(f(x23_signs), x23);
}
