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

// The lp overload used to omit the tanh and exp Jacobian terms of the partial
// correlations and standard deviations; lp must equal log |det J|.
TEST(mathMixMatFun, cov_matrix_constrain_lkj_lp_is_log_det_jacobian) {
  Eigen::VectorXd x4(10);
  x4 << 0.3, -1.2, 0.8, 2.1, -0.4, 1.5, 0.2, -0.7, 1.1, -0.3;
  corr_transform_test::expect_lp_is_log_det_jacobian(
      [](auto&& x) { return stan::math::cov_matrix_constrain_lkj(x, 4); },
      [](auto&& x, auto& lp) {
        stan::math::cov_matrix_constrain_lkj(x, 4, lp);
      },
      x4, false, 1e-12);
  for (int K : {2, 3, 5, 7}) {
    Eigen::VectorXd x((K * (K - 1)) / 2 + K);
    for (int i = 0; i < x.size(); ++i) {
      x(i) = 1.5 * std::sin(0.9 * i + 0.3);
    }
    corr_transform_test::expect_lp_is_log_det_jacobian(
        [K](auto&& x) { return stan::math::cov_matrix_constrain_lkj(x, K); },
        [K](auto&& x, auto& lp) {
          stan::math::cov_matrix_constrain_lkj(x, K, lp);
        },
        x, false, 1e-12);
  }
}
