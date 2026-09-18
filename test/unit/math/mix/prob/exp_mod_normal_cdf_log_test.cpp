#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_exp_mod_normal_lcdf) {
  auto f = [](const auto& p) {
    return stan::math::exp_mod_normal_lcdf(p[0], p[1], p[2], p[3]);
  };

  for (double y : {-40.0, -10.0, -3.0, 0.3, 4.0, 60.0}) {
    stan::test::expect_ad(f, Eigen::Vector4d(y, 0.1, 1.3, 0.7));
  }
  stan::test::expect_ad(
      [](const auto& mu, const auto& sigma, const auto& lambda) {
        return stan::math::exp_mod_normal_lcdf(Eigen::Vector3d(-40, 0.3, 60),
                                               mu, sigma, lambda);
      },
      0.1, 1.3, 0.7);
}

TEST_F(AgradRev, mathMixScalFun_exp_mod_normal_lcdf_tails_and_endpoints) {
  using stan::math::exp_mod_normal_lcdf;
  using stan::math::INFTY;
  using stan::math::NEGATIVE_INFTY;
  // MPFR references for mu = 0, sigma = 1, lambda = 1.
  EXPECT_NEAR(-808.3232158, exp_mod_normal_lcdf(-40.0, 0.0, 1.0, 1.0), 1e-6);
  EXPECT_NEAR(-55.64594635, exp_mod_normal_lcdf(-10.0, 0.0, 1.0, 1.0), 1e-7);
  EXPECT_NEAR(-22.72329719, exp_mod_normal_lcdf(-6.0, 0.0, 1.0, 1.0), 1e-7);
  EXPECT_NEAR(-1.443704555e-26, exp_mod_normal_lcdf(60.0, 0.0, 1.0, 1.0),
              1e-34);
  EXPECT_EQ(NEGATIVE_INFTY, exp_mod_normal_lcdf(NEGATIVE_INFTY, 0.0, 1.0, 1.0));
  EXPECT_EQ(0, exp_mod_normal_lcdf(INFTY, 0.0, 1.0, 1.0));
}
