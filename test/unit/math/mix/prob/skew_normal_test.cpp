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

TEST_F(AgradRev, mathMixScalFun_skew_normal_lpdf_tail) {
  EXPECT_NEAR(-2505.0571524920647,
              stan::math::skew_normal_lpdf(-50.0, 0.0, 1.0, 1.0), 4e-12);
}
