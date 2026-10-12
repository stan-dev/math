#include <test/unit/math/test_ad.hpp>

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

namespace cholesky_corr_constrain_test {
// the strictly lower triangular elements of the Cholesky factor
struct strict_lower_of_constrained {
  int K;
  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& y) const {
    auto L = stan::math::cholesky_corr_constrain(y, K);
    Eigen::Matrix<T, Eigen::Dynamic, 1> v((K * (K - 1)) / 2);
    int k = 0;
    for (int j = 0; j < K; ++j) {
      for (int i = j + 1; i < K; ++i) {
        v(k++) = L(i, j);
      }
    }
    return v;
  }
};
}  // namespace cholesky_corr_constrain_test

// lp is log |det J| of the map from y to the strictly lower triangular
// elements of the Cholesky factor.
TEST(mathMixMatFun, cholesky_corr_constrain_lp_is_log_det_jacobian) {
  Eigen::VectorXd y(6);
  y << 0.3, -1.2, 0.8, 2.1, -0.4, 1.5;
  Eigen::VectorXd fx;
  Eigen::MatrixXd J;
  stan::math::jacobian(
      cholesky_corr_constrain_test::strict_lower_of_constrained{4}, y, fx, J);
  double lp = 0;
  stan::math::cholesky_corr_constrain(y, 4, lp);
  EXPECT_NEAR(std::log(std::fabs(J.determinant())), lp, 1e-12);
}

// lp = sum_n w_n log(1 - tanh(y_n)^2) with weights (1, 1.5, 1) for K = 3,
// and its gradient, at large y against mpmath. Forming 1 - tanh(y)^2 and
// 1 - sum of squares loses digits from |y| = 8.
TEST(mathMixMatFun, cholesky_corr_constrain_lp_large_argument) {
  using stan::math::fvar;
  using stan::math::var;
  Eigen::VectorXd y(3);
  y << 15.0, 12.0, -0.5;
  const double lp_ref = -62.77449311123027;
  Eigen::VectorXd grad_ref(3);
  grad_ref << -1.9999999999996257, -2.9999999997734919, 0.92423431452001952;

  double lp_d = 0;
  stan::math::cholesky_corr_constrain(y, 3, lp_d);
  EXPECT_NEAR(lp_ref, lp_d, 1e-13);

  Eigen::Matrix<var, Eigen::Dynamic, 1> y_v = stan::math::to_var(y);
  var lp_v = 0;
  stan::math::cholesky_corr_constrain(y_v, 3, lp_v);
  lp_v.grad();
  EXPECT_NEAR(lp_ref, lp_v.val(), 1e-13);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(grad_ref(i), y_v(i).adj(), 1e-14);
  }
  stan::math::recover_memory();

  stan::math::var_value<Eigen::VectorXd> y_m(y);
  var lp_m = 0;
  stan::math::cholesky_corr_constrain(y_m, 3, lp_m);
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
    stan::math::cholesky_corr_constrain(y_f, 3, lp_f);
    EXPECT_NEAR(lp_ref, lp_f.val_, 1e-13);
    EXPECT_NEAR(grad_ref(n), lp_f.d_, 1e-14);
  }
}
