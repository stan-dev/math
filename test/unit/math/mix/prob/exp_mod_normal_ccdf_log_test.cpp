#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_exp_mod_normal_lccdf) {
  auto f = [](const auto& p) {
    return stan::math::exp_mod_normal_lccdf(p[0], p[1], p[2], p[3]);
  };

  for (double y : {-40.0, -10.0, -3.0, 0.3, 4.0, 60.0}) {
    stan::test::expect_ad(f, Eigen::Vector4d(y, 0.1, 1.3, 0.7));
  }
  stan::test::expect_ad(
      [](const auto& mu, const auto& sigma, const auto& lambda) {
        return stan::math::exp_mod_normal_lccdf(Eigen::Vector3d(-40, 0.3, 60),
                                                mu, sigma, lambda);
      },
      0.1, 1.3, 0.7);
}

TEST_F(AgradRev, mathMixScalFun_exp_mod_normal_lccdf_tails_and_endpoints) {
  using stan::math::exp_mod_normal_lccdf;
  using stan::math::INFTY;
  using stan::math::NEGATIVE_INFTY;
  // MPFR references for mu = 0, sigma = 1, lambda = 1.
  EXPECT_NEAR(-1.353310396e-10, exp_mod_normal_lccdf(-6.0, 0.0, 1.0, 1.0),
              1e-18);
  EXPECT_FLOAT_EQ(-799.5, exp_mod_normal_lccdf(800.0, 0.0, 1.0, 1.0));
  EXPECT_EQ(0, exp_mod_normal_lccdf(NEGATIVE_INFTY, 0.0, 1.0, 1.0));
  EXPECT_EQ(NEGATIVE_INFTY, exp_mod_normal_lccdf(INFTY, 0.0, 1.0, 1.0));
}
