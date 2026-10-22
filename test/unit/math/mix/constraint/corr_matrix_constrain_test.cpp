#include <test/unit/math/test_ad.hpp>
#include <cmath>

namespace corr_matrix_constrain_test {
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

template <typename T>
auto g1(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  return stan::math::corr_matrix_constrain<false>(x, inv_size(x), lp);
}
template <typename T>
auto g2(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  auto a = stan::math::corr_matrix_constrain<true>(x, inv_size(x), lp);
  return a;
}
template <typename T>
auto g3(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  stan::math::corr_matrix_constrain<true>(x, inv_size(x), lp);
  return lp;
}

template <typename T>
inline void expect_corr_matrix_transform(const T& x) {
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  auto f1 = [](const auto& x) { return g1(x); };
  auto f2 = [](const auto& x) { return g2(x); };
  auto f3 = [](const auto& x) { return g3(x); };
  stan::test::expect_ad(f1, x);
  stan::test::expect_ad(f2, x);
  stan::test::expect_ad(tols, f3, x);
}
}  // namespace corr_matrix_constrain_test

TEST(MathMixMatFun, corr_matrixTransform) {
  // sizes must be (n choose 2)

  Eigen::VectorXd v0(0);
  corr_matrix_constrain_test::expect_corr_matrix_transform(v0);

  Eigen::VectorXd v1(1);
  v1 << -1.7;
  corr_matrix_constrain_test::expect_corr_matrix_transform(v1);

  Eigen::VectorXd v3(3);
  v3 << -1.7, 2.9, 0.01;
  corr_matrix_constrain_test::expect_corr_matrix_transform(v3);

  Eigen::VectorXd v6(6);
  v6 << 1, 2, -3, 1.5, 0.2, 2;
  corr_matrix_constrain_test::expect_corr_matrix_transform(v6);
}

TEST(mathMixMatFun, corr_matrix_constrain) {
  auto f = [](int K) {
    return [K](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      return stan::math::corr_matrix_constrain<false>(x1, K, lp);
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

TEST(mathMixMatFun, corr_matrix_constrain_lp) {
  auto f1 = [](int K) {
    return [K](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      return stan::math::corr_matrix_constrain<true>(x1, K, lp);
    };
  };

  auto f2 = [](int K) {
    return [K](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      stan::math::corr_matrix_constrain<true>(x1, K, lp);
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

namespace corr_matrix_constrain_test {
// the strictly lower triangular elements of the correlation matrix
struct strict_lower_of_constrained {
  int K;
  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& y) const {
    auto Omega = stan::math::corr_matrix_constrain(y, K);
    Eigen::Matrix<T, Eigen::Dynamic, 1> v((K * (K - 1)) / 2);
    int k = 0;
    for (int j = 0; j < K; ++j) {
      for (int i = j + 1; i < K; ++i) {
        v(k++) = Omega(i, j);
      }
    }
    return v;
  }
};
}  // namespace corr_matrix_constrain_test

// lp is log |det J| of the map from y to the strictly lower triangular
// elements of the correlation matrix.
TEST(mathMixMatFun, corr_matrix_constrain_lp_is_log_det_jacobian) {
  Eigen::VectorXd y(6);
  y << 0.3, -1.2, 0.8, 2.1, -0.4, 1.5;
  Eigen::VectorXd fx;
  Eigen::MatrixXd J;
  stan::math::jacobian(
      corr_matrix_constrain_test::strict_lower_of_constrained{4}, y, fx, J);
  double lp = 0;
  stan::math::corr_matrix_constrain(y, 4, lp);
  EXPECT_NEAR(std::log(std::fabs(J.determinant())), lp, 1e-12);
}

// lp = sum_n w_n log(1 - tanh(y_n)^2) with weights (1.5, 1.5, 1) for K = 3,
// and its gradient, at large y against mpmath. Forming 1 - tanh(y)^2 loses
// digits from |y| = 8.
TEST(mathMixMatFun, corr_matrix_constrain_lp_large_argument) {
  using stan::math::fvar;
  using stan::math::var;
  Eigen::VectorXd y(3);
  y << 15.0, 12.0, -0.5;
  const double lp_ref = -77.081345930670418;
  Eigen::VectorXd grad_ref(3);
  grad_ref << -2.9999999999994385, -2.9999999997734919, 0.92423431452001952;

  double lp_d = 0;
  stan::math::corr_matrix_constrain(y, 3, lp_d);
  EXPECT_NEAR(lp_ref, lp_d, 1e-13);

  Eigen::Matrix<var, Eigen::Dynamic, 1> y_v = stan::math::to_var(y);
  var lp_v = 0;
  stan::math::corr_matrix_constrain(y_v, 3, lp_v);
  lp_v.grad();
  EXPECT_NEAR(lp_ref, lp_v.val(), 1e-13);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(grad_ref(i), y_v(i).adj(), 1e-14);
  }
  stan::math::recover_memory();

  stan::math::var_value<Eigen::VectorXd> y_m(y);
  var lp_m = 0;
  stan::math::corr_matrix_constrain(y_m, 3, lp_m);
  lp_m.grad();
  EXPECT_NEAR(lp_ref, lp_m.val(), 1e-13);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(grad_ref(i), y_m.adj()(i), 1e-14);
  }
  stan::math::recover_memory();

  for (int n = 0; n < 3; ++n) {
    Eigen::Matrix<fvar<double>, Eigen::Dynamic, 1> y_f(3);
    for (int i = 0; i < 3; ++i) {
      y_f(i) = fvar<double>(y(i), i == n ? 1.0 : 0.0);
    }
    fvar<double> lp_f = 0;
    stan::math::corr_matrix_constrain(y_f, 3, lp_f);
    EXPECT_NEAR(lp_ref, lp_f.val_, 1e-13);
    EXPECT_NEAR(grad_ref(n), lp_f.d_, 1e-14);
  }
}

namespace corr_matrix_constrain_test {
// strictly lower triangle, row by row
template <typename M>
Eigen::VectorXd strict_lower(const M& Omega) {
  const int K = Omega.rows();
  Eigen::VectorXd v((K * (K - 1)) / 2);
  int k = 0;
  for (int i = 0; i < K; ++i) {
    for (int j = 0; j < i; ++j) {
      v(k++) = stan::math::value_of_rec(Omega(i, j));
    }
  }
  return v;
}
}  // namespace corr_matrix_constrain_test

// Correlations at large |y| against mpmath. With sqrt(1 - z^2) from a
// rounded z, the small elements of the Cholesky factor had relative errors
// up to 3.8e-4, which reached the correlations close to 1.
TEST(mathMixMatFun, corr_matrix_constrain_values_large_argument) {
  using corr_matrix_constrain_test::strict_lower;
  using stan::math::fvar;
  using stan::math::var;
  Eigen::VectorXd y(6);
  y << 15.0, 12.0, -0.5, 9.0, 2.0, -14.0;
  Eigen::VectorXd ref(6);
  ref << 0.99999999999981285, 0.99999999992449731, 0.99999999993182827,
      -0.46211715726000976, -0.46211663421720721, -0.46210665234659286;
  auto expect_ref = [&ref](const Eigen::VectorXd& v, const char* mode) {
    for (int i = 0; i < ref.size(); ++i) {
      EXPECT_NEAR(ref(i), v(i), 4e-16) << mode << ", element " << i;
    }
  };
  expect_ref(strict_lower(stan::math::corr_matrix_constrain(y, 4)), "double");
  Eigen::Matrix<var, Eigen::Dynamic, 1> y_v = stan::math::to_var(y);
  expect_ref(strict_lower(stan::math::corr_matrix_constrain(y_v, 4)), "var");
  stan::math::var_value<Eigen::VectorXd> y_m(y);
  expect_ref(
      strict_lower(stan::math::corr_matrix_constrain(y_m, 4).val().eval()),
      "var_value");
  Eigen::Matrix<fvar<double>, Eigen::Dynamic, 1> y_f(6);
  for (int i = 0; i < 6; ++i) {
    y_f(i) = fvar<double>(y(i), 1.0);
  }
  expect_ref(strict_lower(stan::math::corr_matrix_constrain(y_f, 4)), "fvar");
  stan::math::recover_memory();
}

// K = 40 with unconstrained values in [-1.9, 1.9], as default inits: the
// diagonal must be 1, lp finite, and the var_value gradient must equal the
// Eigen<var> gradient and be finite.
TEST(mathMixMatFun, corr_matrix_constrain_long_rows) {
  using stan::math::var;
  const int K = 40;
  Eigen::VectorXd y((K * (K - 1)) / 2);
  for (int i = 0; i < y.size(); ++i) {
    y(i) = 1.9 * std::sin(1.3 * i + 0.7);
  }
  Eigen::MatrixXd W(K, K);
  for (int i = 0; i < W.size(); ++i) {
    W(i) = std::cos(0.7 * i);
  }
  double lp = 0;
  Eigen::MatrixXd Omega = stan::math::corr_matrix_constrain(y, K, lp);
  EXPECT_LT((Omega.diagonal().array() - 1).abs().maxCoeff(), 1e-14);
  EXPECT_TRUE(std::isfinite(lp));

  Eigen::Matrix<var, Eigen::Dynamic, 1> y_v = stan::math::to_var(y);
  var lp_v = 0;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> Omega_v
      = stan::math::corr_matrix_constrain(y_v, K, lp_v);
  var f_v = (W.array() * Omega_v.array()).sum() + lp_v;
  f_v.grad();
  Eigen::VectorXd g_v = stan::math::value_of(y_v.adj());
  stan::math::recover_memory();

  stan::math::var_value<Eigen::VectorXd> y_m(y);
  var lp_m = 0;
  stan::math::var_value<Eigen::MatrixXd> Omega_m
      = stan::math::corr_matrix_constrain(y_m, K, lp_m);
  var f_m = stan::math::sum(stan::math::elt_multiply(W, Omega_m)) + lp_m;
  f_m.grad();
  Eigen::VectorXd g_m = y_m.adj();
  stan::math::recover_memory();

  EXPECT_TRUE(g_v.allFinite());
  EXPECT_TRUE(g_m.allFinite());
  EXPECT_LT((g_v - g_m).cwiseAbs().maxCoeff(),
            1e-12 * g_v.cwiseAbs().maxCoeff());
}
