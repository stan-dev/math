#include <test/unit/math/test_ad.hpp>

namespace cov_matrix_constrain_test {
// easier than fiddling the quadratic equation
template <typename T>
int inv_size(const T& x) {
  int dof = x.size();
  for (int k = 0; k < 10; ++k) {
    if (dof == k + (k * (k - 1)) / 2) {
      return k;
    }
  }
  return -1;
}

template <typename T>
auto g1(const T& x) {
  stan::scalar_type_t<T> lp = 0.0;
  return stan::math::cov_matrix_constrain<false>(x, inv_size(x), lp);
}
template <typename T>
auto g2(const T& x) {
  typename stan::scalar_type<T>::type lp = 0;
  auto a = stan::math::cov_matrix_constrain<true>(x, inv_size(x), lp);
  return a;
}
template <typename T>
auto g3(const T& x) {
  typename stan::scalar_type<T>::type lp = 0;
  stan::math::cov_matrix_constrain<true>(x, inv_size(x), lp);
  return lp;
}

template <typename T>
void expect_cov_matrix_transform(const T& x) {
  using stan::test::relative_tolerance;
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = relative_tolerance(1e-3, 1e-3);
  tols.hessian_fvar_hessian_ = relative_tolerance(1e-3, 1e-3);

  auto f1 = [](const auto& x) { return g1(x); };
  auto f2 = [](const auto& x) { return g2(x); };
  auto f3 = [](const auto& x) { return g3(x); };
  stan::test::expect_ad(f1, x);
  stan::test::expect_ad(f2, x);
  stan::test::expect_ad(tols, f3, x);
}
}  // namespace cov_matrix_constrain_test

TEST(MathMixMatFun, cov_matrixTransform) {
  // sizes must be n + (n choose 2)

  Eigen::VectorXd v0(0);
  cov_matrix_constrain_test::expect_cov_matrix_transform(v0);

  // 1 x 1
  Eigen::VectorXd v1(1);
  v1 << -1.7;
  cov_matrix_constrain_test::expect_cov_matrix_transform(v1);

  // 2 x 2
  Eigen::VectorXd v3(3);
  v3 << -1.7, 2.9, 0.01;
  cov_matrix_constrain_test::expect_cov_matrix_transform(v3);

  // 3 x 3
  Eigen::VectorXd v6(6);
  v6 << 1, 2, -3, 1.5, 0.2, 2;
  cov_matrix_constrain_test::expect_cov_matrix_transform(v6);

  // 4 x 4
  Eigen::VectorXd v10(10);
  v10 << 1, 2, -3, 1.7, 9.8, -12.2, 0.4, 0.2, 1.2, 2.7;
  cov_matrix_constrain_test::expect_cov_matrix_transform(v10);
}

TEST(mathMixMatFun, cov_matrix_constrain) {
  auto f = [](int K) {
    return [K](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      return stan::math::cov_matrix_constrain<false>(x1, K, lp);
    };
  };

  Eigen::VectorXd x1(10);
  x1 << -0.9, 0.2, 0.99, 0.1, 0.2, 0.3, -0.1, -0.2, -0.3, -0.4;
  Eigen::VectorXd x2(6);
  x2 << -0.3, 0.2, -0.99, 0.1, 0.2, 0.3;
  stan::test::expect_ad(f(4), x1);
  stan::test::expect_ad(f(3), x2);
  stan::test::expect_ad_matvar(f(4), x1);
  stan::test::expect_ad_matvar(f(3), x2);
}

TEST(mathMixMatFun, cov_matrix_constrain_lp) {
  auto f1 = [](int K) {
    return [K](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      return stan::math::cov_matrix_constrain<true>(x1, K, lp);
    };
  };

  auto f2 = [](int K) {
    return [K](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      stan::math::cov_matrix_constrain<true>(x1, K, lp);
      return lp;
    };
  };

  Eigen::VectorXd x1(10);
  x1 << -0.9, 0.0, 0.99, 0.1, 0.2, 0.3, -0.1, -0.2, -0.3, -0.4;
  Eigen::VectorXd x2(6);
  x2 << -0.3, 0.2, -0.99, 0.1, 0.2, 0.3;
  stan::test::expect_ad(f1(4), x1);
  stan::test::expect_ad(f1(3), x2);

  stan::test::expect_ad_matvar(f1(4), x1);
  stan::test::expect_ad_matvar(f1(3), x2);

  stan::test::expect_ad(f2(4), x1);
  stan::test::expect_ad(f2(3), x2);

  stan::test::expect_ad_matvar(f2(4), x1);
  stan::test::expect_ad_matvar(f2(3), x2);
}
