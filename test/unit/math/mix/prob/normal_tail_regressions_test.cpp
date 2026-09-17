#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <cmath>

TEST_F(AgradRev, normal_lcdf_infinite_endpoints) {
  using namespace stan::math;
  EXPECT_EQ(0, normal_lcdf(INFTY, 0.0, 1.0));
  EXPECT_EQ(NEGATIVE_INFTY, normal_lcdf(NEGATIVE_INFTY, 0.0, 1.0));
  var sigma = 1;
  auto lp = normal_lcdf(INFTY, 0.0, sigma);
  lp.grad();
  EXPECT_EQ(0, sigma.adj());
}

TEST_F(AgradRev, lognormal_shared_normal_tail) {
  using namespace stan::math;
  EXPECT_NEAR(-1254.8313611394199, lognormal_lcdf(std::exp(-50.0), 0.0, 1.0),
              2e-12);
  stan::test::expect_ad([](const auto& y, const auto& m,
                           const auto& s) { return lognormal_lcdf(y, m, s); },
                        1.0, 50.0, 1.0);
  const Eigen::Vector3d y(1, std::exp(-50.0), 2);
  stan::test::expect_ad(
      [&](const auto& m, const auto& s) { return lognormal_lcdf(y, m, s); },
      0.0, 1.0);
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

TEST_F(AgradRev, skew_normal_shared_normal_tail) {
  using namespace stan::math;
  EXPECT_NEAR(-2505.0571524920647, skew_normal_lpdf(-50.0, 0.0, 1.0, 1.0),
              4e-12);
  const auto f
      = [](const auto& p) { return skew_normal_lpdf(p[0], p[1], p[2], p[3]); };
  stan::test::expect_ad(f, Eigen::Vector4d(-50, 0, 1, 1));
  stan::test::expect_ad(f, Eigen::Vector4d(0.3, 0.1, 1.2, 0.5));
  stan::test::expect_ad(
      [](const auto& m, const auto& s, const auto& a) {
        return skew_normal_lpdf(Eigen::Vector3d(-50, 0, 2), m, s, a);
      },
      0.1, 1.2, 0.5);
}

TEST_F(AgradRev, lognormal_lccdf_cdf_shared_normal_tail) {
  using namespace stan::math;
  EXPECT_NEAR(-1254.8313611394199, lognormal_lccdf(std::exp(50.0), 0.0, 1.0),
              2e-12);
  EXPECT_NEAR(1,
              lognormal_cdf(std::exp(-20.0), 0.0, 1.0) / 2.7536241186062337e-89,
              1e-12);
  EXPECT_EQ(0, lognormal_lccdf(0.0, 0.0, 1.0));
  EXPECT_EQ(0, lognormal_cdf(0.0, 0.0, 1.0));
  const Eigen::Vector3d y(1, std::exp(-50.0), std::exp(50.0));
  stan::test::expect_ad(
      [&](const auto& m, const auto& s) { return lognormal_lccdf(y, m, s); },
      0.0, 1.0);
  stan::test::expect_ad(
      [&](const auto& m, const auto& s) { return lognormal_cdf(y, m, s); }, 0.0,
      1.0);
  stan::test::expect_ad([](const auto& y, const auto& m,
                           const auto& s) { return lognormal_lccdf(y, m, s); },
                        2.0, 0.5, 1.5);
  stan::test::expect_ad([](const auto& y, const auto& m,
                           const auto& s) { return lognormal_cdf(y, m, s); },
                        2.0, 0.5, 1.5);
}

TEST_F(AgradRev, exp_mod_normal_log_space) {
  using namespace stan::math;
  // MPFR references for mu = 0, sigma = 1, lambda = 1.
  EXPECT_NEAR(-808.3232158, exp_mod_normal_lcdf(-40.0, 0.0, 1.0, 1.0), 1e-6);
  EXPECT_NEAR(-55.64594635, exp_mod_normal_lcdf(-10.0, 0.0, 1.0, 1.0), 1e-7);
  EXPECT_NEAR(-22.72329719, exp_mod_normal_lcdf(-6.0, 0.0, 1.0, 1.0), 1e-7);
  EXPECT_NEAR(-1.353310396e-10, exp_mod_normal_lccdf(-6.0, 0.0, 1.0, 1.0),
              1e-18);
  EXPECT_NEAR(-1.443704555e-26, exp_mod_normal_lcdf(60.0, 0.0, 1.0, 1.0),
              1e-34);
  EXPECT_FLOAT_EQ(-799.5, exp_mod_normal_lccdf(800.0, 0.0, 1.0, 1.0));
  EXPECT_EQ(NEGATIVE_INFTY, exp_mod_normal_lcdf(NEGATIVE_INFTY, 0.0, 1.0, 1.0));
  EXPECT_EQ(0, exp_mod_normal_lcdf(INFTY, 0.0, 1.0, 1.0));
  EXPECT_EQ(0, exp_mod_normal_lccdf(NEGATIVE_INFTY, 0.0, 1.0, 1.0));
  EXPECT_EQ(NEGATIVE_INFTY, exp_mod_normal_lccdf(INFTY, 0.0, 1.0, 1.0));
  EXPECT_EQ(0, exp_mod_normal_cdf(NEGATIVE_INFTY, 0.0, 1.0, 1.0));
  EXPECT_EQ(1, exp_mod_normal_cdf(INFTY, 0.0, 1.0, 1.0));
  EXPECT_NEAR(-1.6127097401997972, exp_mod_normal_lpdf(0.0, 0.0, 2.0, 20.0),
              1e-12);

  const auto lcdf = [](const auto& p) {
    return exp_mod_normal_lcdf(p[0], p[1], p[2], p[3]);
  };
  const auto lccdf = [](const auto& p) {
    return exp_mod_normal_lccdf(p[0], p[1], p[2], p[3]);
  };
  const auto cdf = [](const auto& p) {
    return exp_mod_normal_cdf(p[0], p[1], p[2], p[3]);
  };
  const auto lpdf = [](const auto& p) {
    return exp_mod_normal_lpdf(p[0], p[1], p[2], p[3]);
  };
  for (double y : {-40.0, -10.0, -3.0, 0.3, 4.0, 60.0}) {
    const Eigen::Vector4d p(y, 0.1, 1.3, 0.7);
    stan::test::expect_ad(lcdf, p);
    stan::test::expect_ad(lccdf, p);
    stan::test::expect_ad(cdf, p);
    stan::test::expect_ad(lpdf, p);
  }
  stan::test::expect_ad(lpdf, Eigen::Vector4d(0, 0, 2, 20));
  const Eigen::Vector3d y(-40, 0.3, 60);
  stan::test::expect_ad(
      [&](const auto& m, const auto& s, const auto& l) {
        return exp_mod_normal_lcdf(y, m, s, l);
      },
      0.1, 1.3, 0.7);
  stan::test::expect_ad(
      [&](const auto& m, const auto& s, const auto& l) {
        return exp_mod_normal_lccdf(y, m, s, l);
      },
      0.1, 1.3, 0.7);
  stan::test::expect_ad(
      [&](const auto& m, const auto& s, const auto& l) {
        return exp_mod_normal_lpdf(y, m, s, l);
      },
      0.1, 1.3, 0.7);
}
