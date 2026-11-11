#include <test/unit/math/mix/constraint/corr_transform_helpers.hpp>
#include <cmath>

namespace cholesky_corr_constrain_test {
// easier than fiddling the quadratic equation
template <typename T>
int inv_size(const T& x) {
  int k_choose_2 = x.size();
  for (int k = 0; k < 10; ++k) {
    if (k_choose_2 == (k * (k - 1)) / 2) {
      return k;
    }
  }
  return -1;
}

// K (K - 1) / 2 deterministic values in [-1.9, 1.9], like default inits
inline Eigen::VectorXd long_rows_input(int K) {
  Eigen::VectorXd y((K * (K - 1)) / 2);
  for (int i = 0; i < y.size(); ++i) {
    y(i) = 1.9 * std::sin(1.3 * i + 0.7);
  }
  return y;
}

template <typename T>
auto g1(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  return stan::math::cholesky_corr_constrain<false>(x, inv_size(x), lp);
}
template <typename T>
auto g2(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  auto a = stan::math::cholesky_corr_constrain<true>(x, inv_size(x), lp);
  return a;
}
template <typename T>
auto g3(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  stan::math::cholesky_corr_constrain<true>(x, inv_size(x), lp);
  return lp;
}

template <typename T>
inline void expect_cholesky_corr_transform(const T& x) {
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  auto f1 = [](const auto& x) { return cholesky_corr_constrain_test::g1(x); };
  auto f2 = [](const auto& x) { return cholesky_corr_constrain_test::g2(x); };
  auto f3 = [](const auto& x) { return cholesky_corr_constrain_test::g3(x); };
  stan::test::expect_ad(f1, x);
  stan::test::expect_ad(f2, x);
  stan::test::expect_ad(tols, f3, x);
}
}  // namespace cholesky_corr_constrain_test

TEST(MathMixMatFun, cholesky_corrTransform) {
  // sizes must be (n choose 2)

  Eigen::VectorXd v0(0);
  cholesky_corr_constrain_test::expect_cholesky_corr_transform(v0);

  Eigen::VectorXd v1(1);
  v1 << -1.7;
  cholesky_corr_constrain_test::expect_cholesky_corr_transform(v1);

  Eigen::VectorXd v3(3);
  v3 << -1.7, 2.9, 0.01;
  cholesky_corr_constrain_test::expect_cholesky_corr_transform(v3);

  Eigen::VectorXd v6(6);
  v6 << 1, 2, -3, 1.5, 0.2, 2;
  cholesky_corr_constrain_test::expect_cholesky_corr_transform(v6);
}

TEST(mathMixMatFun, cholesky_corr_constrain) {
  auto f = [](int K) {
    return [K](const auto& x1) {
      return stan::math::cholesky_corr_constrain(x1, K);
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

TEST(mathMixMatFun, cholesky_corr_constrain_lp) {
  auto f1 = [](int K) {
    return [K](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      return stan::math::cholesky_corr_constrain(x1, K, lp);
    };
  };

  auto f2 = [](int K) {
    return [K](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      stan::math::cholesky_corr_constrain(x1, K, lp);
      return lp;
    };
  };

  Eigen::VectorXd x1(6);
  x1 << -0.9, 0.0, 0.99, 0.1, 0.2, 0.3;
  Eigen::VectorXd x2(3);
  x2 << -0.3, 0.2, -0.99;
  stan::test::expect_ad(f1(4), x1);
  stan::test::expect_ad(f1(3), x2);

  stan::test::expect_ad_matvar(f1(4), x1);
  stan::test::expect_ad_matvar(f1(3), x2);

  stan::test::expect_ad(f2(4), x1);
  stan::test::expect_ad(f2(3), x2);

  stan::test::expect_ad_matvar(f2(4), x1);
  stan::test::expect_ad_matvar(f2(3), x2);
}

// The Jacobian is a closed-form weighted sum of log_cosh(y); the weights
// are only right if lp equals log |det J| of the autodiff Jacobian.
TEST(mathMixMatFun, cholesky_corr_constrain_lp_is_log_det_jacobian) {
  Eigen::VectorXd y4(6);
  y4 << 0.3, -1.2, 0.8, 2.1, -0.4, 1.5;
  corr_transform_test::expect_lp_is_log_det_jacobian(
      [](auto&& y) { return stan::math::cholesky_corr_constrain(y, 4); },
      [](auto&& y, auto& lp) { stan::math::cholesky_corr_constrain(y, 4, lp); },
      y4, true, 1e-12);
  for (int K : {2, 3, 5, 7}) {
    Eigen::VectorXd y((K * (K - 1)) / 2);
    for (int i = 0; i < y.size(); ++i) {
      y(i) = 1.5 * std::sin(0.9 * i + 0.3);
    }
    corr_transform_test::expect_lp_is_log_det_jacobian(
        [K](auto&& y) { return stan::math::cholesky_corr_constrain(y, K); },
        [K](auto&& y, auto& lp) {
          stan::math::cholesky_corr_constrain(y, K, lp);
        },
        y, true, 1e-12);
  }
}

// The old Jacobian terms log1m(square(z)) and log1m(sum_sqs) cancel at
// large |y|; the log_cosh form must stay accurate there in every mode.
TEST(mathMixMatFun, cholesky_corr_constrain_lp_large_argument) {
  Eigen::VectorXd y(3);
  y << 15.0, 12.0, -0.5;
  Eigen::VectorXd grad_ref(3);
  grad_ref << -1.9999999999996257, -2.9999999997734919, 0.92423431452001952;
  corr_transform_test::expect_lp_all_modes(
      [](auto&& y, auto& lp) { stan::math::cholesky_corr_constrain(y, 3, lp); },
      y, -62.77449311123027, grad_ref, 1e-13, 1e-14);
}

// The old sqrt(1 - sum_sqs) recursion loses the small elements of L when
// partial correlations are near 1; the tanh / log_cosh factor must not.
TEST(mathMixMatFun, cholesky_corr_constrain_values_large_argument) {
  Eigen::VectorXd y(6);
  y << 15.0, 12.0, -0.5, 9.0, 2.0, -14.0;
  Eigen::VectorXd ref(10);
  ref << 1.0, 0.99999999999981285, 6.1180464100359433e-7, 0.99999999992449731,
      -5.6786918924293556e-6, 1.0897607083695929e-5, 0.99999996954004097,
      0.00023794090595878951, -6.5605200973176475e-5, 1.091052174635906e-10;
  corr_transform_test::expect_values_all_modes(
      [](auto&& y) { return stan::math::cholesky_corr_constrain(y, 4); }, y,
      false, ref, 0.0, 4e-15);
}

// #3218: with default inits at K = 40 a diagonal element of L rounded to 0,
// lkj_corr_cholesky_lpdf rejected it, and initialization failed.
TEST(mathMixMatFun, cholesky_corr_constrain_long_rows) {
  const int K = 40;
  Eigen::VectorXd y = cholesky_corr_constrain_test::long_rows_input(K);
  double lp = 0;
  Eigen::MatrixXd L = stan::math::cholesky_corr_constrain(y, K, lp);
  EXPECT_GT(L.diagonal().minCoeff(), 0.0);
  EXPECT_LT((L.rowwise().squaredNorm().array() - 1).abs().maxCoeff(), 1e-14);
  EXPECT_TRUE(std::isfinite(lp));
  EXPECT_TRUE(std::isfinite(stan::math::lkj_corr_cholesky_lpdf(L, 1.0)));
}

// The var_value overload has a hand-written reverse pass; Eigen<var>
// differentiates the forward formula and is the oracle for it at K = 40.
TEST(mathMixMatFun, cholesky_corr_constrain_long_rows_var_value_gradient) {
  using stan::math::var;
  const int K = 40;
  Eigen::VectorXd y = cholesky_corr_constrain_test::long_rows_input(K);

  Eigen::Matrix<var, Eigen::Dynamic, 1> y_v = stan::math::to_var(y);
  var lp_v = 0;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> L_v
      = stan::math::cholesky_corr_constrain(y_v, K, lp_v);
  var f_v = 2.0 * stan::math::sum(stan::math::log(L_v.diagonal())) + lp_v;
  f_v.grad();
  Eigen::VectorXd g_v = stan::math::value_of(y_v.adj());
  stan::math::recover_memory();

  stan::math::var_value<Eigen::VectorXd> y_m(y);
  var lp_m = 0;
  stan::math::var_value<Eigen::MatrixXd> L_m
      = stan::math::cholesky_corr_constrain(y_m, K, lp_m);
  var f_m = 2.0 * stan::math::sum(stan::math::log(L_m.diagonal())) + lp_m;
  f_m.grad();
  Eigen::VectorXd g_m = y_m.adj();
  stan::math::recover_memory();

  EXPECT_TRUE(g_v.allFinite());
  EXPECT_TRUE(g_m.allFinite());
  EXPECT_LT((g_v - g_m).cwiseAbs().maxCoeff(),
            1e-12 * g_v.cwiseAbs().maxCoeff());
}
