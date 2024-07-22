#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(mathMixScalFun, logSumExp) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::log_sum_exp(x1, x2);
  };
  stan::test::expect_ad(f, 0.0, 0.0);
  stan::test::expect_ad(f, 0.5, -1.0);
  stan::test::expect_ad(f, 0.5, 0.0);
  stan::test::expect_ad(f, 0.5, 1.2);
  stan::test::expect_ad(f, 0.5, 1.4);
  stan::test::expect_ad(f, 1.0, 2.0);
  stan::test::expect_ad(f, 1.4, 0.5);
  stan::test::expect_ad(f, 1.2, 0.6);
  stan::test::expect_ad(f, 2.0, 1.0);
  stan::test::expect_ad(f, 3.0, 5.0);
  stan::test::expect_ad(f, 3.0, 6.0);
  stan::test::expect_ad(f, 3.4, 0.9);
  stan::test::expect_ad(f, 5.0, 2.0);

  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 0.01;
  tols.gradient_fvar_grad_ = 0.01;
  tols.hessian_grad_ = 0.01;
  tols.hessian_hessian_ = 3.0;
  tols.hessian_fvar_grad_ = 0.01;
  tols.hessian_fvar_hessian_ = 3.0;
  tols.grad_hessian_hessian_ = 3.0;
  tols.grad_hessian_grad_hessian_ = 3.0;
  stan::test::expect_ad(tols, f, -2.0, 15.0);
  stan::test::expect_ad(tols, f, 2.0, -15.0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);

  stan::test::expect_value(f, 1000.0, 10.0);
}

TEST(MathMixMatFun, logSumExp) {
  auto f = [](const auto& x) { return stan::math::log_sum_exp(x); };

  Eigen::VectorXd x0(0);
  stan::test::expect_ad(f, x0);
  stan::test::expect_ad_matvar(f, x0);

  Eigen::VectorXd x1(1);
  x1 << 0;

  Eigen::VectorXd x2(2);
  x2 << 5, 2;

  Eigen::VectorXd x2b(2);
  x2b << 4.9, -std::numeric_limits<double>::infinity();

  Eigen::VectorXd x2c(2);
  x2c << std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity();

  Eigen::VectorXd x4(4);
  x4 << 1, 2, 3, 4;

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  for (const auto& x : std::vector<Eigen::VectorXd>{x1, x2, x2b, x2c, x4}) {
    stan::test::expect_ad(tols, f, x);
    stan::test::expect_ad_matvar(tols, f, x);
    Eigen::RowVectorXd rx = x;
    stan::test::expect_ad(tols, f, rx);
    stan::test::expect_ad_matvar(tols, f, rx);
    std::vector<double> stx
        = std::vector<double>(x.data(), x.data() + x.size());
    stan::test::expect_ad(tols, f, stx);
  }

  Eigen::MatrixXd x23(2, 2);
  x23 << 1, 2, 3, 4;
  stan::test::expect_ad(f, x23);
  stan::test::expect_ad_matvar(f, x23);

  std::vector<Eigen::VectorXd> stvx{x2, x2b, x2c};
  stan::test::expect_ad(tols, f, stvx);
  stan::test::expect_ad_matvar(tols, f, stvx);
  std::vector<Eigen::RowVectorXd> strx{x2.transpose(), x2b.transpose(),
                                       x2c.transpose()};
  stan::test::expect_ad(tols, f, strx);
  stan::test::expect_ad_matvar(tols, f, strx);
  std::vector<std::vector<double>> ststx{
      std::vector<double>(x2.data(), x2.data() + x2.size()),
      std::vector<double>(x2b.data(), x2b.data() + x2b.size()),
      std::vector<double>(x2c.data(), x2c.data() + x2c.size())};
  stan::test::expect_ad(tols, f, ststx);
}

TEST(mathMixScalFun, logSumExp_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::log_sum_exp;
    return log_sum_exp(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}
