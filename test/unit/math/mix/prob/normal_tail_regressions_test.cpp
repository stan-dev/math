#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <cmath>

TEST_F(AgradRev, normal_standardization_overflow) {
  using namespace stan::math;
  var y = -1e308, mu = 1e308, sigma = 1e308;
  auto lp = normal_lcdf(y, mu, sigma);
  lp.grad();
  const auto ref = internal::std_normal_lcdf_value_grad<true>(-2.0);
  EXPECT_NEAR(ref.first, lp.val(), 1e-14);
  EXPECT_NEAR(ref.second / 1e308, y.adj(), 1e-320);
  EXPECT_NEAR(-ref.second / 1e308, mu.adj(), 1e-320);
  EXPECT_NEAR(2 * (ref.second / 1e308), sigma.adj(), 1e-320);
  EXPECT_NEAR(std_normal_cdf(-2.0), normal_cdf(-1e308, 1e308, 1e308), 1e-15);
  EXPECT_NEAR(std_normal_lpdf(-2.0) - std::log(1e308),
              normal_lpdf(-1e308, 1e308, 1e308), 1e-12);
  EXPECT_TRUE(std::isfinite(normal_lpdf(-1.5e154, 0.0, 1.0)));
  EXPECT_DOUBLE_EQ(normal_lpdf(3.0, 1.0, 4.0), normal_lpdf(3, 1, 4));
  const double near_mu = std::nextafter(1e308, 0.0);
  EXPECT_EQ((1e308 - near_mu) / 1e292,
            internal::normal_standardize(1e308, near_mu, 1e292));
  const Eigen::Vector3d ys(-1e308, 1e308, 1.1);
  const Eigen::Vector3d ms(1e308, 1e308, 1.0);
  const Eigen::Vector3d ss(1e308, 1e308, 0.1);
  EXPECT_NEAR(std_normal_lcdf(Eigen::Vector3d(-2, 0, 1)),
              normal_lcdf(ys, ms, ss), 1e-14);
  // Infinite observations retain their endpoint behavior.
  EXPECT_EQ(0, normal_lcdf(INFTY, 0.0, 1.0));
  EXPECT_EQ(NEGATIVE_INFTY, normal_lcdf(NEGATIVE_INFTY, 0.0, 1.0));
}

TEST_F(AgradRev, normal_lpdf_extreme_scale_derivatives) {
  using namespace stan::math;
  var sigma = 1e154;
  auto lp = normal_lpdf(1.5e308, 0.0, sigma);
  lp.grad();
  EXPECT_TRUE(std::isfinite(lp.val()));
  EXPECT_NEAR(2.25e154, sigma.adj(), 1e140);
  set_zero_all_adjoints();
  var y = 0;
  auto centered = normal_lpdf(y, 0.0, 1e-310);
  centered.grad();
  EXPECT_TRUE(std::isfinite(centered.val()));
  EXPECT_EQ(0, y.adj());
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
