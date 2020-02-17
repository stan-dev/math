#include <stan/math/mix.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <gtest/gtest.h>
#include <vector>

template <typename F>
void expect_match_autodiff(const F& f, Eigen::VectorXd x) {
  double fx_fd;
  Eigen::VectorXd grad_fx_fd;
  stan::math::finite_diff_gradient_auto(f, x, fx_fd, grad_fx_fd);

  double fx;
  Eigen::VectorXd grad_fx;
  stan::math::gradient(f, x, fx, grad_fx);

  EXPECT_FLOAT_EQ(fx, fx_fd);
  EXPECT_EQ(grad_fx.size(), grad_fx_fd.size());
  stan::test::expect_near_rel("expect_match_autodiff", grad_fx, grad_fx_fd,
                              stan::test::relative_tolerance(1e-7, 1e-9));
}

TEST(MathMixMatFunctor, FiniteDiffGradientAuto) {
  auto norm_fun
      = [](const auto& x) { return stan::math::normal_lpdf(x(0), x(1), x(2)); };

  for (double frac = 1e-40; frac <= 1e+40; frac *= 1e10) {
    Eigen::VectorXd y(3);
    y << 1 * frac, 2 * frac, 3;
    expect_match_autodiff(norm_fun, y);
  }

  for (double frac = 1; frac <= 1e+10; frac *= 10) {
    Eigen::VectorXd y(3);
    y << 0, 0, frac;
    expect_match_autodiff(norm_fun, y);
  }

  auto log_fun
      = [](const auto& x) { return stan::math::sum(stan::math::log(x)); };
  Eigen::VectorXd y(0);
  expect_match_autodiff(log_fun, y);

  Eigen::VectorXd z(1);
  z << 2;
  expect_match_autodiff(log_fun, z);

  Eigen::VectorXd w(5);
  w << 1, 2, 3, 4, 5;
  expect_match_autodiff(log_fun, w);
}
