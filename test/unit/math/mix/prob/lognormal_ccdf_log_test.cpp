#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <cmath>

TEST_F(AgradRev, mathMixScalFun_lognormal_lccdf) {
  auto f = [](const auto& y, const auto& mu, const auto& sigma) {
    return stan::math::lognormal_lccdf(y, mu, sigma);
  };

  stan::test::expect_ad(f, 2.0, 0.5, 1.5);
  const Eigen::Vector3d y(1, std::exp(-50.0), std::exp(50.0));
  stan::test::expect_ad(
      [&](const auto& mu, const auto& sigma) {
        return stan::math::lognormal_lccdf(y, mu, sigma);
      },
      0.0, 1.0);
}

TEST_F(AgradRev, mathMixScalFun_lognormal_lccdf_tail_and_endpoints) {
  using stan::math::lognormal_lccdf;
  EXPECT_NEAR(-1254.8313611394199, lognormal_lccdf(std::exp(50.0), 0.0, 1.0),
              2e-12);
  EXPECT_EQ(0, lognormal_lccdf(0.0, 0.0, 1.0));
}
