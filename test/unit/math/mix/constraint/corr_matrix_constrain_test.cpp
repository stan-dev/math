#include <test/unit/math/mix/constraint/corr_transform_helpers.hpp>
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

// The Jacobian is a closed-form weighted sum of log_cosh(y) in column order;
// the weights are only right if lp equals log |det J| of the autodiff Jacobian.
TEST(mathMixMatFun, corr_matrix_constrain_lp_is_log_det_jacobian) {
  Eigen::VectorXd y4(6);
  y4 << 0.3, -1.2, 0.8, 2.1, -0.4, 1.5;
  corr_transform_test::expect_lp_is_log_det_jacobian(
      [](auto&& y) { return stan::math::corr_matrix_constrain(y, 4); },
      [](auto&& y, auto& lp) { stan::math::corr_matrix_constrain(y, 4, lp); },
      y4, true, 1e-12);
  for (int K : {2, 3, 5, 7}) {
    Eigen::VectorXd y((K * (K - 1)) / 2);
    for (int i = 0; i < y.size(); ++i) {
      y(i) = 1.5 * std::sin(0.9 * i + 0.3);
    }
    corr_transform_test::expect_lp_is_log_det_jacobian(
        [K](auto&& y) { return stan::math::corr_matrix_constrain(y, K); },
        [K](auto&& y, auto& lp) {
          stan::math::corr_matrix_constrain(y, K, lp);
        },
        y, true, 1e-12);
  }
}

// The old read_corr_L Jacobian terms cancel at large |y|; the log_cosh form
// must stay accurate there in every mode.
TEST(mathMixMatFun, corr_matrix_constrain_lp_large_argument) {
  Eigen::VectorXd y(3);
  y << 15.0, 12.0, -0.5;
  Eigen::VectorXd grad_ref(3);
  grad_ref << -2.9999999999994385, -2.9999999997734919, 0.92423431452001952;
  corr_transform_test::expect_lp_all_modes(
      [](auto&& y, auto& lp) { stan::math::corr_matrix_constrain(y, 3, lp); },
      y, -77.081345930670418, grad_ref, 1e-13, 1e-14);
}

// The old sqrt(1 - sum_sqs) recursion put ~1e-10 errors in the correlations
// when the partial correlations are near 1; the log_cosh factor must not.
TEST(mathMixMatFun, corr_matrix_constrain_values_large_argument) {
  Eigen::VectorXd y(6);
  y << 15.0, 12.0, -0.5, 9.0, 2.0, -14.0;
  Eigen::VectorXd ref(6);
  ref << 0.99999999999981285, 0.99999999992449731, 0.99999999993182827,
      -0.46211715726000976, -0.46211663421720721, -0.46210665234659286;
  corr_transform_test::expect_values_all_modes(
      [](auto&& y) { return stan::math::corr_matrix_constrain(y, 4); }, y, true,
      ref, 4e-16, 0.0);
}
