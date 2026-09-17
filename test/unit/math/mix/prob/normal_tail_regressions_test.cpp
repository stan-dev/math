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
