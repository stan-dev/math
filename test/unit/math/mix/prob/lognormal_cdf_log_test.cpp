#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <cmath>

TEST_F(AgradRev, mathMixScalFun_lognormal_lcdf) {
  auto f = [](const auto& y, const auto& mu, const auto& sigma) {
    return stan::math::lognormal_lcdf(y, mu, sigma);
  };

  stan::test::expect_ad(f, 2.0, 0.5, 1.5);
  stan::test::expect_ad(f, 1.0, 50.0, 1.0);
  const Eigen::Vector3d y(1, std::exp(-50.0), 2);
  stan::test::expect_ad(
      [&](const auto& mu, const auto& sigma) {
        return stan::math::lognormal_lcdf(y, mu, sigma);
      },
      0.0, 1.0);
}

TEST_F(AgradRev, mathMixScalFun_lognormal_lcdf_tail_and_endpoints) {
  using stan::math::INFTY;
  using stan::math::lognormal_lcdf;
  using stan::math::NEGATIVE_INFTY;
  using stan::math::var;
  EXPECT_NEAR(-1254.8313611394199, lognormal_lcdf(std::exp(-50.0), 0.0, 1.0),
              2e-12);
  EXPECT_EQ(NEGATIVE_INFTY, lognormal_lcdf(0.0, 0.0, 1.0));
  var sigma = 1;
  auto lp = lognormal_lcdf(INFTY, 0.0, sigma);
  lp.grad();
  EXPECT_EQ(0, lp.val());
  EXPECT_EQ(0, sigma.adj());
  EXPECT_THROW(
      lognormal_lcdf(Eigen::VectorXd::Ones(2), Eigen::VectorXd::Ones(3), 1),
      std::invalid_argument);
}
