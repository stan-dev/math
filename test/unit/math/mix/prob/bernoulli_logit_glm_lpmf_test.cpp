#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <limits>
#include <vector>

TEST_F(AgradRev, mathMixScalFun_bernoulli_logit_glm_lpmf) {
  auto f = [](const auto y) {
    return [=](const auto& x, const auto& alpha, const auto& beta) {
      return stan::math::bernoulli_logit_glm_lpmf(y, x, alpha, beta);
    };
  };

  std::vector<int> y{0, 1};
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);
  Eigen::RowVectorXd x_rowvec = x.row(0);
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(2);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(2);

  stan::test::expect_ad(f(y[0]), x, alpha, beta);
  stan::test::expect_ad(f(y[0]), x, alpha[0], beta);
  stan::test::expect_ad(f(y[1]), x_rowvec, alpha, beta);
  stan::test::expect_ad(f(y[1]), x_rowvec, alpha[0], beta);
  stan::test::expect_ad(f(y), x, alpha, beta);
  stan::test::expect_ad(f(y), x, alpha[0], beta);
  stan::test::expect_ad(f(y), x_rowvec, alpha, beta);
  stan::test::expect_ad(f(y), x_rowvec, alpha[0], beta);
}

namespace bernoulli_logit_glm_test {
// theta = alpha + x beta = J v with v = (alpha, beta), J = [I, x]. With
// s = 2y - 1, t = s theta, p = inv_logit(t), q = 1 - p, the n-th summand is
// log(p_n) with theta_n derivatives s q, -p q and -s p q (q - p); the
// derivatives in v follow by the chain rule through J. Each tolerance is
// relative to the same sum taken over absolute values.
template <typename T_y, typename T_x>
void expect_derivatives(const T_y& y, const T_x& x,
                        const Eigen::VectorXd& alpha,
                        const Eigen::VectorXd& beta) {
  using stan::math::fvar;
  using stan::math::var;
  const int N = alpha.size();
  const int K = beta.size();
  const int M = N + K;
  Eigen::MatrixXd J(N, M);
  J << Eigen::MatrixXd::Identity(N, N), x.replicate(N / x.rows(), 1);
  Eigen::VectorXd v(M);
  v << alpha, beta;
  Eigen::VectorXd theta = J * v;
  Eigen::ArrayXd lp(N), d1(N), d2(N), d3(N), w(N);
  for (int n = 0; n < N; ++n) {
    double s = 2 * stan::get(y, n) - 1, t = s * theta(n);
    double p = stan::math::inv_logit(t), q = stan::math::inv_logit(-t);
    lp(n) = stan::math::log_inv_logit(t);
    d1(n) = s * q;
    d2(n) = -p * q;
    d3(n) = -s * p * q * (q - p);
    w(n) = p * q;
  }
  const double tiny = std::numeric_limits<double>::min();
  Eigen::MatrixXd A = J.cwiseAbs();
  Eigen::VectorXd grad_exact = J.transpose() * d1.matrix();
  Eigen::VectorXd grad_tol
      = (1e-13 * A.transpose() * d1.abs().matrix()).array() + tiny;
  Eigen::MatrixXd H_exact = J.transpose() * d2.matrix().asDiagonal() * J;
  Eigen::MatrixXd H_tol
      = (1e-13 * A.transpose() * w.matrix().asDiagonal() * A).array() + tiny;

  auto f = [&](const auto& u) {
    return stan::math::bernoulli_logit_glm_lpmf(y, x, u.head(N), u.tail(K));
  };
  double fx;
  Eigen::VectorXd grad, grad_fwd;
  Eigen::MatrixXd H, H_fwd, H_ffv;
  std::vector<Eigen::MatrixXd> grad_H;
  stan::math::hessian(f, v, fx, grad, H);
  stan::math::hessian<double>(f, v, fx, grad_fwd, H_fwd);  // fvar<fvar<double>>
  stan::math::grad_hessian(f, v, fx, H_ffv, grad_H);
  EXPECT_NEAR(fx, lp.sum(), 1e-13 * lp.abs().sum() + tiny);
  for (int i = 0; i < M; ++i) {
    EXPECT_NEAR(grad(i), grad_exact(i), grad_tol(i)) << i;
    EXPECT_NEAR(grad_fwd(i), grad_exact(i), grad_tol(i)) << i;
    for (int j = 0; j < M; ++j) {
      EXPECT_NEAR(H(i, j), H_exact(i, j), H_tol(i, j)) << i << ", " << j;
      EXPECT_NEAR(H_fwd(i, j), H_exact(i, j), H_tol(i, j)) << i << ", " << j;
      EXPECT_NEAR(H_ffv(i, j), H_exact(i, j), H_tol(i, j)) << i << ", " << j;
      Eigen::ArrayXd Jij = J.col(i).array() * J.col(j).array();
      for (int k = 0; k < M; ++k) {
        Eigen::ArrayXd Jijk = Jij * J.col(k).array();
        EXPECT_NEAR(grad_H[i](j, k), (d3 * Jijk).sum(),
                    1e-13 * (w * Jijk.abs()).sum() + tiny)
            << i << ", " << j << ", " << k;
      }
    }
  }

  // derivative of the value itself
  stan::math::nested_rev_autodiff nested;
  Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1> v_fv(M);
  for (int k = 0; k < M; ++k) {
    v_fv(k) = fvar<var>(v(k), 0);
  }
  f(v_fv).val_.grad();
  for (int i = 0; i < M; ++i) {
    EXPECT_NEAR(v_fv(i).val_.adj(), grad_exact(i), grad_tol(i)) << i;
    Eigen::Matrix<fvar<fvar<double>>, Eigen::Dynamic, 1> v_ff(M);
    for (int k = 0; k < M; ++k) {
      v_ff(k) = fvar<fvar<double>>(fvar<double>(v(k), i == k), 0);
    }
    EXPECT_NEAR(f(v_ff).val_.d_, grad_exact(i), grad_tol(i)) << i;
  }
}
}  // namespace bernoulli_logit_glm_test

// t = (2y - 1) theta in the tails (exp(-t) overflows at t < -709, curvature
// was 0 below -20) and at t = +-0, where the value's own derivative
// must not switch branch. Matrix and row-vector x, vector and scalar y.
TEST_F(AgradRev, mathMixScalFun_bernoulli_logit_glm_lpmf_tails) {
  using bernoulli_logit_glm_test::expect_derivatives;
  std::vector<int> y{1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0};
  Eigen::MatrixXd x(13, 2);
  x << 0.5, -1, 1, 0.25, -0.5, 2, 1.5, -0.75, 0.25, 0.125, -1, -0.25, 0.75, 0.5,
      -0.625, 1.25, 0.375, -1, 1.125, 0.625, -0.5, -1.25, 0.75, 0.5, -0.25, 1;
  Eigen::VectorXd beta(2);  // dyadic x and beta: theta = 0 exactly
  beta << 0.5, -0.25;
  Eigen::VectorXd theta(13);
  theta << -2500, 720, -30, 20.5, -0.25, 0, 0, -0.0, -0.0, 21, 21, 40, -19.5;
  Eigen::VectorXd alpha = theta - x * beta;
  expect_derivatives(y, x, alpha, beta);
  expect_derivatives(1, x, alpha, beta);

  Eigen::RowVectorXd x_row = x.row(0);
  alpha = theta.array() - x_row.dot(beta);
  expect_derivatives(y, x_row, alpha, beta);
}
