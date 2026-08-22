#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST_F(AgradRev, mathMixScalFun_inv_gaussian_lpdf) {
  auto f = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lpdf(y, mu, lambda);
  };

  stan::test::expect_ad(f, 0.3, 0.5, 2.0);
  stan::test::expect_ad(f, 1.2, 0.5, 2.0);
  stan::test::expect_ad(f, 0.3, 1.0, 5.0);
  stan::test::expect_ad(f, 5.0, 1.0, 0.5);
  stan::test::expect_ad(f, 1e-3, 1.0, 2.0);
  // out of support and invalid parameters
  stan::test::expect_ad(f, -1.0, 1.0, 2.0);
  stan::test::expect_ad(f, 1.0, 0.0, 2.0);
  stan::test::expect_ad(f, 1.0, 1.0, 0.0);
}

TEST_F(AgradRev, mathMixScalFun_inv_gaussian_cdf) {
  auto f = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_cdf(y, mu, lambda);
  };

  stan::test::expect_ad(f, 0.3, 0.5, 2.0);
  stan::test::expect_ad(f, 1.2, 0.5, 2.0);
  stan::test::expect_ad(f, 0.3, 1.0, 5.0);
  stan::test::expect_ad(f, 5.0, 1.0, 0.5);
  // out of support and invalid parameters
  stan::test::expect_ad(f, -1.0, 1.0, 2.0);
  stan::test::expect_ad(f, 1.0, 0.0, 2.0);
}

TEST_F(AgradRev, mathMixScalFun_inv_gaussian_lcdf) {
  auto f = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lcdf(y, mu, lambda);
  };

  stan::test::expect_ad(f, 0.3, 0.5, 2.0);
  stan::test::expect_ad(f, 1.2, 0.5, 2.0);
  stan::test::expect_ad(f, 0.3, 1.0, 5.0);
  stan::test::expect_ad(f, 5.0, 1.0, 0.5);
  // out of support and invalid parameters
  stan::test::expect_ad(f, -1.0, 1.0, 2.0);
  stan::test::expect_ad(f, 1.0, 0.0, 2.0);
}

TEST_F(AgradRev, mathMixScalFun_inv_gaussian_lccdf) {
  auto f = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lccdf(y, mu, lambda);
  };

  stan::test::expect_ad(f, 0.3, 0.5, 2.0);
  stan::test::expect_ad(f, 1.2, 0.5, 2.0);
  stan::test::expect_ad(f, 0.3, 1.0, 5.0);
  stan::test::expect_ad(f, 5.0, 1.0, 0.5);
  // out of support and invalid parameters
  stan::test::expect_ad(f, -1.0, 1.0, 2.0);
  stan::test::expect_ad(f, 1.0, 0.0, 2.0);
}

// crosses the internal log_Phi branch at z = -30
TEST_F(AgradRev, mathMixScalFun_inv_gaussian_lcdf_tails) {
  auto f = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lcdf(y, mu, lambda);
  };

  stan::test::expect_ad(f, 0.1, 0.1, 100.0);
  stan::test::expect_ad(f, 0.05, 1.0, 100.0);
  stan::test::expect_ad(f, 0.15, 0.1, 10.0);
}

TEST_F(AgradRev, mathMixScalFun_inv_gaussian_lccdf_tails) {
  auto f = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lccdf(y, mu, lambda);
  };

  stan::test::expect_ad(f, 0.1, 0.1, 100.0);
  stan::test::expect_ad(f, 0.15, 0.1, 10.0);
  stan::test::expect_ad(f, 0.5, 0.1, 50.0);
}

TEST_F(AgradRev, mathMixScalFun_inv_gaussian_lpdf_vectorized) {
  Eigen::VectorXd y(3);
  y << 0.3, 1.2, 4.0;
  Eigen::VectorXd mu(3);
  mu << 0.5, 1.0, 2.0;
  Eigen::VectorXd lambda(3);
  lambda << 2.0, 5.0, 0.7;

  auto f = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lpdf(y, mu, lambda);
  };
  stan::test::expect_ad(f, y, mu, lambda);
}

TEST_F(AgradRev, mathMixScalFun_inv_gaussian_lccdf_vectorized) {
  Eigen::VectorXd y(3);
  y << 0.3, 1.2, 4.0;
  Eigen::VectorXd mu(3);
  mu << 0.5, 1.0, 2.0;
  Eigen::VectorXd lambda(3);
  lambda << 2.0, 5.0, 0.7;

  auto f = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lccdf(y, mu, lambda);
  };
  stan::test::expect_ad(f, y, mu, lambda);
}
