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
void expect_cholesky_corr_transform(const T& x) {
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
