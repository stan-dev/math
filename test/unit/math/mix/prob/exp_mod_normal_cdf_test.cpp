#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_exp_mod_normal_cdf) {
  auto f = [](const auto& p) {
    return stan::math::exp_mod_normal_cdf(p[0], p[1], p[2], p[3]);
  };

  for (double y : {-40.0, -10.0, -3.0, 0.3, 4.0, 60.0}) {
    stan::test::expect_ad(f, Eigen::Vector4d(y, 0.1, 1.3, 0.7));
  }
}

TEST_F(AgradRev, mathMixScalFun_exp_mod_normal_cdf_endpoints) {
  using stan::math::exp_mod_normal_cdf;
  using stan::math::INFTY;
  using stan::math::NEGATIVE_INFTY;
  EXPECT_EQ(0, exp_mod_normal_cdf(NEGATIVE_INFTY, 0.0, 1.0, 1.0));
  EXPECT_EQ(1, exp_mod_normal_cdf(INFTY, 0.0, 1.0, 1.0));
}
