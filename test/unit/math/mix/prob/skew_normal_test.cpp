#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_skew_normal_lpdf) {
  auto f = [](const auto& p) {
    return stan::math::skew_normal_lpdf(p[0], p[1], p[2], p[3]);
  };

  stan::test::expect_ad(f, Eigen::Vector4d(0.3, 0.1, 1.2, 0.5));
  stan::test::expect_ad(f, Eigen::Vector4d(-50, 0, 1, 1));
  stan::test::expect_ad(
      [](const auto& mu, const auto& sigma, const auto& alpha) {
        return stan::math::skew_normal_lpdf(Eigen::Vector3d(-50, 0, 2), mu,
                                            sigma, alpha);
      },
      0.1, 1.2, 0.5);
}

TEST_F(AgradRev, mathMixScalFun_skew_normal_lpdf_vector_tails) {
  // For alpha = +/-1, alpha z reaches the deep lower tail, crosses
  // +/-4 sqrt(2) and passes 40.
  Eigen::VectorXd y(6);
  y << -1e3, -50, -5.7, 0, 5.7, 50;
  for (double alpha : {-1.0, 1.0}) {
    stan::test::expect_ad(
        [&y](const auto& mu, const auto& sigma, const auto& alpha) {
          return stan::math::skew_normal_lpdf(y, mu, sigma, alpha);
        },
        0.0, 1.0, alpha);
  }
  Eigen::VectorXd alpha(7);
  alpha << -300, -30, -1, 0, 1, 30, 300;
  stan::test::expect_ad(
      [](const auto& y, const auto& sigma, const auto& alpha) {
        return stan::math::skew_normal_lpdf(y, 0.1, sigma, alpha);
      },
      0.3, 1.2, alpha);
  Eigen::VectorXd mu(4);
  mu << -10, 0, 0.3, 8;
  Eigen::VectorXd sigma(4);
  sigma << 0.5, 1, 2, 3;
  stan::test::expect_ad(
      [](const auto& mu, const auto& sigma, const auto& alpha) {
        return stan::math::skew_normal_lpdf(0.3, mu, sigma, alpha);
      },
      mu, sigma, 4.0);
}

TEST_F(AgradRev, mathMixScalFun_skew_normal_lpdf_tail) {
  EXPECT_NEAR(-2505.0571524920647,
              stan::math::skew_normal_lpdf(-50.0, 0.0, 1.0, 1.0), 4e-12);
}
