#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <string>
#include <vector>

template <typename F>
void test_grad_hessian_finite_diff(const std::string& msg, const F& f,
                                   Eigen::VectorXd& x) {
  double fx;
  Eigen::VectorXd grad_fx;
  Eigen::MatrixXd hess_fx;
  std::vector<Eigen::MatrixXd> grad_hess_fx;
  stan::math::finite_diff_grad_hessian_auto(f, x, fx, hess_fx, grad_hess_fx);

  double fx_ad;
  Eigen::VectorXd grad_fx_ad;
  Eigen::MatrixXd hess_fx_ad;
  std::vector<Eigen::MatrixXd> grad_hess_fx_ad;
  stan::math::grad_hessian(f, x, fx_ad, hess_fx_ad, grad_hess_fx_ad);

  EXPECT_FLOAT_EQ(fx_ad, fx) << msg;

  EXPECT_EQ(hess_fx_ad.rows(), hess_fx.rows()) << msg;
  EXPECT_EQ(hess_fx_ad.cols(), hess_fx.cols()) << msg;
  for (int i = 0; i < hess_fx_ad.size(); ++i)
    EXPECT_NEAR(hess_fx_ad(i), hess_fx(i), 1e-4) << msg;

  EXPECT_EQ(grad_hess_fx_ad.size(), grad_hess_fx.size());
  for (size_t i = 0; i < grad_hess_fx_ad.size(); ++i) {
    EXPECT_EQ(grad_hess_fx_ad[i].rows(), grad_hess_fx[i].rows());
    EXPECT_EQ(grad_hess_fx_ad[i].cols(), grad_hess_fx[i].cols());
    for (int j = 0; j < grad_hess_fx_ad[i].size(); ++j)
      EXPECT_NEAR(grad_hess_fx_ad[i](j), grad_hess_fx[i](j), 1e-3);
  }
}

TEST(MixMatFunctor, FiniteDiffGradHessianAuto) {
  auto norm_fun
      = [](const auto& x) { return stan::math::normal_lpdf(x(0), x(1), x(2)); };
  Eigen::VectorXd x(3);
  for (const auto& scale : std::vector<double>{1, 1e10, 1e20, 1e30}) {
    x << 1 * scale, 2 * scale, 1.5 * scale;
    test_grad_hessian_finite_diff("norm_fun({1, 2, 3})", norm_fun, x);
  }
  for (const auto& scale : std::vector<double>{1e-10, 1e-20, 1e-30}) {
    x << scale, -scale, 1;  // finite diff fails with small sigma
    test_grad_hessian_finite_diff("norm_fun({1, 2, 3})", norm_fun, x);
  }

  auto log_fun
      = [](const auto& x) { return stan::math::sum(stan::math::log(x)); };
  Eigen::VectorXd y(0);
  test_grad_hessian_finite_diff("log_fun({})", log_fun, y);

  Eigen::VectorXd z(1);
  z << 2;
  test_grad_hessian_finite_diff("log_fun({2})", log_fun, z);

  Eigen::VectorXd w(5);
  w << 1, 2, 3, 4, 5;
  test_grad_hessian_finite_diff("log_fun({1, 2, 3, 4, 5})", log_fun, w);
}
