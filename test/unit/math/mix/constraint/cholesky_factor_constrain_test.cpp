#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, cholesky_factor_constrain) {
  auto f = [](int M, int N) {
    return [M, N](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      return stan::math::cholesky_factor_constrain<false>(x1, M, N, lp);
    };
  };

  Eigen::VectorXd x1(14);
  x1 << -0.9, 0.2, 0.99, 0.1, 0.2, 0.3, -0.1, -0.2, -0.3, -0.4, -0.4, -0.5,
      -0.6, -0.7;
  Eigen::VectorXd x2(12);
  x2 << -0.3, 0.2, -0.99, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, -0.9;
  stan::test::expect_ad(f(5, 4), x1);
  stan::test::expect_ad(f(5, 3), x2);
  stan::test::expect_ad_matvar(f(5, 4), x1);
  stan::test::expect_ad_matvar(f(5, 3), x2);
}

TEST(mathMixMatFun, cholesky_factor_constrain_lp) {
  auto f1 = [](int M, int N) {
    return [M, N](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      return stan::math::cholesky_factor_constrain<true>(x1, M, N, lp);
    };
  };

  auto f2 = [](int M, int N) {
    return [M, N](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      stan::math::cholesky_factor_constrain<true>(x1, M, N, lp);
      return lp;
    };
  };

  Eigen::VectorXd x1(14);
  x1 << -0.9, 0.2, 0.99, 0.1, 0.2, 0.3, -0.1, -0.2, -0.3, -0.4, -0.4, -0.5,
      -0.6, -0.7;
  Eigen::VectorXd x2(12);
  x2 << -0.3, 0.2, -0.99, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, -0.9;
  stan::test::expect_ad(f1(5, 4), x1);
  stan::test::expect_ad(f1(5, 3), x2);

  stan::test::expect_ad_matvar(f1(5, 4), x1);
  stan::test::expect_ad_matvar(f1(5, 3), x2);

  stan::test::expect_ad(f2(5, 4), x1);
  stan::test::expect_ad(f2(5, 3), x2);

  stan::test::expect_ad_matvar(f2(5, 4), x1);
  stan::test::expect_ad_matvar(f2(5, 3), x2);
}
