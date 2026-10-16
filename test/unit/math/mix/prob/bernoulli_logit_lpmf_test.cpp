#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <limits>
#include <vector>

namespace bernoulli_logit_lpmf_test {
// With s = 2y - 1, t = s theta, p = inv_logit(t), q = 1 - p, the log density
// is log(p) and its theta derivatives are s q, -p q and -s p q (q - p).
struct exact {
  double lp, d1, d2, d3, scale;
  exact(int y, double theta) {
    double s = 2 * y - 1, t = s * theta;
    double p = stan::math::inv_logit(t), q = stan::math::inv_logit(-t);
    lp = stan::math::log_inv_logit(t);
    d1 = s * q;
    d2 = -p * q;
    d3 = -s * p * q * (q - p);
    scale = p * q;
  }
};

inline double tol(double x) {
  return 1e-13 * std::fabs(x) + std::numeric_limits<double>::min();
}
}  // namespace bernoulli_logit_lpmf_test

// t = (2y - 1) theta in the tails (exp(-t) overflows at t < -709, curvature
// was 0 below -20) and at t = +-0, where the value's own derivative
// must not switch branch. Closed forms in every autodiff mode, scalar and
// vector theta.
TEST_F(AgradRev, mathMixScalFun_bernoulli_logit_lpmf_tails) {
  using bernoulli_logit_lpmf_test::exact;
  using bernoulli_logit_lpmf_test::tol;
  using stan::math::fvar;
  using stan::math::var;
  std::vector<int> y{1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0};
  Eigen::VectorXd theta(13);
  theta << -2500, 720, -30, 20.5, -0.25, 0, 0, -0.0, -0.0, 21, 21, 40, -19.5;
  const int N = theta.size();

  for (int n = 0; n < N; ++n) {
    exact e(y[n], theta(n));
    auto f = [&](const auto& v) {
      return stan::math::bernoulli_logit_lpmf(y[n], v(0));
    };
    Eigen::VectorXd v = theta.segment(n, 1);
    double fx;
    Eigen::VectorXd grad;
    Eigen::MatrixXd H;
    stan::math::hessian(f, v, fx, grad, H);
    EXPECT_NEAR(fx, e.lp, tol(e.lp)) << n;
    EXPECT_NEAR(grad(0), e.d1, tol(e.d1)) << n;
    EXPECT_NEAR(H(0, 0), e.d2, tol(e.d2)) << n;
    stan::math::hessian<double>(f, v, fx, grad, H);  // fvar<fvar<double>>
    EXPECT_NEAR(grad(0), e.d1, tol(e.d1)) << n;
    EXPECT_NEAR(H(0, 0), e.d2, tol(e.d2)) << n;
    std::vector<Eigen::MatrixXd> grad_H;
    stan::math::grad_hessian(f, v, fx, H, grad_H);
    EXPECT_NEAR(H(0, 0), e.d2, tol(e.d2)) << n;
    EXPECT_NEAR(grad_H[0](0, 0), e.d3, tol(e.scale)) << n;

    // derivative of the value itself
    fvar<fvar<double>> theta_ff(fvar<double>(theta(n), 1), 0);
    EXPECT_NEAR(stan::math::bernoulli_logit_lpmf(y[n], theta_ff).val_.d_, e.d1,
                tol(e.d1))
        << n;
    stan::math::nested_rev_autodiff nested;
    fvar<var> theta_fv(theta(n), 0);
    stan::math::bernoulli_logit_lpmf(y[n], theta_fv).val_.grad();
    EXPECT_NEAR(theta_fv.val_.adj(), e.d1, tol(e.d1)) << n;
  }

  auto f
      = [&](const auto& v) { return stan::math::bernoulli_logit_lpmf(y, v); };
  double fx;
  Eigen::VectorXd grad;
  Eigen::MatrixXd H;
  stan::math::hessian(f, theta, fx, grad, H);
  double lp = 0, lp_tol = 0;
  for (int i = 0; i < N; ++i) {
    exact e(y[i], theta(i));
    lp += e.lp;
    lp_tol += tol(e.lp);
    EXPECT_NEAR(grad(i), e.d1, tol(e.d1)) << i;
    for (int j = 0; j < N; ++j) {
      EXPECT_NEAR(H(i, j), i == j ? e.d2 : 0, tol(e.d2)) << i << ", " << j;
    }
  }
  EXPECT_NEAR(fx, lp, lp_tol);
}
