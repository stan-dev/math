#include <test/unit/math/mix/constraint/corr_transform_helpers.hpp>
#include <cmath>

TEST(mathMixMatFun, cov_matrix_constrain_lkj) {
  auto f = [](int K) {
    return [K](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      return stan::math::cov_matrix_constrain_lkj<false>(x1, K, lp);
    };
  };

  Eigen::VectorXd x1(6);
  x1 << -0.9, 0.2, 0.99, 0.1, 0.2, 0.3;
  Eigen::VectorXd x2(3);
  x2 << -0.3, 0.2, -0.99;
  stan::test::expect_ad(f(4), x1);
  stan::test::expect_ad(f(3), x2);
  stan::test::expect_ad_matvar(f(4), x1);
  stan::test::expect_ad_matvar(f(3), x2);
}

TEST(mathMixMatFun, cov_matrix_constrain_lkj_lp) {
  auto f1 = [](int K) {
    return [K](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      return stan::math::cov_matrix_constrain_lkj<true>(x1, K, lp);
    };
  };

  auto f2 = [](int K) {
    return [K](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      stan::math::cov_matrix_constrain_lkj<true>(x1, K, lp);
      return lp;
    };
  };

  Eigen::VectorXd x1(6);
  x1 << -0.9, 0.0, 0.99, 0.1, 0.2, 0.3;
  Eigen::VectorXd x2(3);
  x2 << -0.3, 0.2, -0.99;
  Eigen::VectorXd x3(1);
  x3 << 0.1;
  stan::test::expect_ad(f1(4), x1);
  stan::test::expect_ad(f1(3), x2);

  stan::test::expect_ad_matvar(f1(4), x1);
  stan::test::expect_ad_matvar(f1(3), x2);

  stan::test::expect_ad(f2(4), x1);
  stan::test::expect_ad(f2(3), x2);

  stan::test::expect_ad_matvar(f2(4), x1);
  stan::test::expect_ad_matvar(f2(3), x2);
}

// lp is log |det J| of the map from x (6 partial correlations, then 4 log
// standard deviations) to the lower triangle of the covariance matrix,
// diagonal included.
TEST(mathMixMatFun, cov_matrix_constrain_lkj_lp_is_log_det_jacobian) {
  Eigen::VectorXd x(10);
  x << 0.3, -1.2, 0.8, 2.1, -0.4, 1.5, 0.2, -0.7, 1.1, -0.3;
  corr_transform_test::expect_lp_is_log_det_jacobian(
      [](const auto& x) { return stan::math::cov_matrix_constrain_lkj(x, 4); },
      [](const auto& x, auto& lp) {
        stan::math::cov_matrix_constrain_lkj(x, 4, lp);
      },
      x, false, 1e-12);
}

// K = 40 with unconstrained values in [-1.9, 1.9], as default inits: the
// covariance matrix must be the correlation matrix scaled by exp(x) with
// the diagonal exp(2 x), lp finite, and the var_value gradient must equal
// the Eigen<var> gradient and be finite. (The matrix has a condition number
// near 1e26, so a floating-point Cholesky factorization of it is not a
// meaningful check.)
TEST(mathMixMatFun, cov_matrix_constrain_lkj_long_rows) {
  using stan::math::var;
  const int K = 40;
  const int k_choose_2 = (K * (K - 1)) / 2;
  Eigen::VectorXd x(k_choose_2 + K);
  for (int i = 0; i < x.size(); ++i) {
    x(i) = 1.9 * std::sin(1.3 * i + 0.7);
  }
  Eigen::MatrixXd W(K, K);
  for (int i = 0; i < W.size(); ++i) {
    W(i) = std::cos(0.7 * i);
  }
  double lp = 0;
  Eigen::MatrixXd Sigma = stan::math::cov_matrix_constrain_lkj(x, K, lp);
  EXPECT_TRUE(std::isfinite(lp));
  EXPECT_LT((Sigma.diagonal().array() / (2 * x.tail(K).array()).exp() - 1)
                .abs()
                .maxCoeff(),
            1e-14);
  Eigen::VectorXd sds = x.tail(K).array().exp();
  Eigen::MatrixXd expected
      = sds.asDiagonal()
        * stan::math::corr_matrix_constrain(x.head(k_choose_2).eval(), K)
        * sds.asDiagonal();
  EXPECT_LT((Sigma - expected).cwiseAbs().maxCoeff(),
            1e-14 * expected.cwiseAbs().maxCoeff());

  Eigen::Matrix<var, Eigen::Dynamic, 1> x_v = stan::math::to_var(x);
  var lp_v = 0;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> Sigma_v
      = stan::math::cov_matrix_constrain_lkj(x_v, K, lp_v);
  var f_v = (W.array() * Sigma_v.array()).sum() + lp_v;
  f_v.grad();
  Eigen::VectorXd g_v = stan::math::value_of(x_v.adj());
  stan::math::recover_memory();

  stan::math::var_value<Eigen::VectorXd> x_m(x);
  var lp_m = 0;
  stan::math::var_value<Eigen::MatrixXd> Sigma_m
      = stan::math::cov_matrix_constrain_lkj(x_m, K, lp_m);
  var f_m = stan::math::sum(stan::math::elt_multiply(W, Sigma_m)) + lp_m;
  f_m.grad();
  Eigen::VectorXd g_m = x_m.adj();
  stan::math::recover_memory();

  EXPECT_TRUE(g_v.allFinite());
  EXPECT_TRUE(g_m.allFinite());
  EXPECT_LT((g_v - g_m).cwiseAbs().maxCoeff(),
            1e-12 * g_v.cwiseAbs().maxCoeff());
}
