#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

// 1 / sigma overflows for sigma below 1 / DBL_MAX (about 5.6e-309). At the
// median the value is log(1/2) for every sigma, and the y and mu partials are
// +-phi(0) / Phi(0) / sigma, which overflows to +-inf for such sigma. The
// sigma partial is 0 where 1 / sigma is finite.

namespace {

template <typename F>
void expect_median(const F& f, double y0, double sigma0, double sign) {
  using stan::math::var;
  var y = y0, mu = 0.0, sigma = sigma0;
  var lp = f(y, mu, sigma);
  lp.grad();
  EXPECT_DOUBLE_EQ(lp.val(), std::log(0.5)) << "sigma = " << sigma0;
  // phi(0) / Phi(0) = sqrt(2 / pi)
  const double slope = std::sqrt(2.0 / stan::math::pi()) / sigma0;
  if (std::isfinite(slope)) {
    EXPECT_NEAR(y.adj(), sign * slope, 1e-14 * slope);
    EXPECT_NEAR(mu.adj(), -sign * slope, 1e-14 * slope);
    EXPECT_EQ(sigma.adj(), 0.0) << "sigma = " << sigma0;
  } else {
    EXPECT_EQ(y.adj(), sign * std::numeric_limits<double>::infinity());
    EXPECT_EQ(mu.adj(), -sign * std::numeric_limits<double>::infinity());
  }
  stan::math::recover_memory();
}

const std::vector<double> sigmas{1.0, 1e-300, 1e-308, 1e-310, 4.9e-324};

}  // namespace

TEST_F(AgradRev, ProbNormalLcdf_median_with_subnormal_sigma) {
  for (double s : sigmas) {
    expect_median(
        [](auto& y, auto& mu, auto& sigma) {
          return stan::math::normal_lcdf(y, mu, sigma);
        },
        0.0, s, 1.0);
    expect_median(
        [](auto& y, auto& mu, auto& sigma) {
          return stan::math::normal_lccdf(y, mu, sigma);
        },
        0.0, s, -1.0);
  }
}

TEST_F(AgradRev, ProbLognormalLcdf_median_with_subnormal_sigma) {
  for (double s : sigmas) {
    expect_median(
        [](auto& y, auto& mu, auto& sigma) {
          return stan::math::lognormal_lcdf(y, mu, sigma);
        },
        1.0, s, 1.0);
    expect_median(
        [](auto& y, auto& mu, auto& sigma) {
          return stan::math::lognormal_lccdf(y, mu, sigma);
        },
        1.0, s, -1.0);
  }
}

TEST_F(AgradRev, ProbNormalLcdf_median_with_subnormal_sigma_vector) {
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> sigma(3);
  sigma << 1.0, 1e-310, 2.0;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y(3);
  y << 0.0, 0.0, 1.0;
  var lp = stan::math::normal_lcdf(y, 0.0, sigma);
  lp.grad();
  const double expected
      = 2 * std::log(0.5) + std::log(0.5 * std::erfc(-0.5 / std::sqrt(2.0)));
  EXPECT_NEAR(lp.val(), expected, 1e-15);
  EXPECT_TRUE(std::isfinite(sigma(2).adj()));
  stan::math::recover_memory();
}
