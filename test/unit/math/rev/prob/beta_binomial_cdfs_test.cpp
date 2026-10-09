#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <vector>

namespace {

// Values and shape gradients of beta_binomial_cdf, beta_binomial_lcdf and
// beta_binomial_lccdf at points where the 3F2 series of these functions
// ends at a zero numerator parameter. Reference values: sums of the pmf and
// of its closed-form shape derivatives (mpmath, 50 digits).
struct beta_binomial_cdfs_point {
  int n;
  int N;
  double alpha;
  double beta;
  // rows: cdf, lcdf, lccdf; columns: value, d/dalpha, d/dbeta
  double ref[3][3];
};

const std::vector<beta_binomial_cdfs_point> beta_binomial_cdfs_points = {
    // large shapes: the 3F2 terms grow before the series ends
    {3,
     10,
     298.30333970,
     146.75306521,
     {{0.019577975830469777, -0.00025213357314530317, 0.0005048815375647927},
      {-3.9333500266686489, -0.012878429074005716, 0.025788239904711231},
      {-0.019772163107297, 0.00025716841006185449, -0.00051496348013239942}}},
    // beta = 1: a denominator parameter of 3F2 equals the numerator
    // parameter that ends the series
    {3,
     10,
     2.5,
     1.0,
     {{0.10178821044348033, -0.078609438975422664, 0.12549085581810137},
      {-2.2848609925400993, -0.77228432087497911, 1.232862384271726},
      {-0.1073493926558798, 0.08751771006505609, -0.13971187784126151}}},
    {8,
     10,
     3.0,
     1.0,
     {{0.57692307692307692, -0.092455621301775148, 0.3774924861463323},
      {-0.55004633691927198, -0.16025641025641026, 0.65432030932030932},
      {-0.8602012652231115, 0.21853146853146853, -0.89225496725496725}}},
    // alpha + beta = 1: the sums of the 3F2 numerator and denominator
    // parameters are equal
    {3,
     10,
     0.5,
     0.5,
     {{0.4080352783203125, -0.61252655150398375, 0.49626754881843688},
      {-0.89640164184463105, -1.5011607673372379, 1.2162368677074559},
      {-0.52430823763107531, 1.0347348905623168, -0.83833973654762416}}},
    // N - n = 1: a 3F2 numerator parameter is zero
    {0,
     1,
     100.0,
     0.5,
     {{0.0049751243781094527, -0.000049503725155317938, 0.0099007450310635875},
      {-5.3033049080590758, -0.0099502487562189055, 1.9900497512437811},
      {-0.0049875415110390736, 0.000049751243781094527,
       -0.0099502487562189055}}},
    // the terms of the 117-term 3F2 series fall below 1e-6 before the end
    // of the series; the second point has the mirrored series
    {0,
     117,
     0.1,
     10.0,
     {{0.77231348114687892, -1.9916580324442808, 0.0074694843556493395},
      {-0.25836474771713782, -2.5788207522762463, 0.0096715705966406531},
      {-1.4797855134046491, 8.7473691568409674, -0.032806001836533192}}},
    {116,
     117,
     10.0,
     0.1,
     {{0.22768651885312108, -0.0074694843556493395, 1.9916580324442808},
      {-1.4797855134046491, -0.032806001836533192, 8.7473691568409674},
      {-0.25836474771713782, 0.0096715705966406531, -2.5788207522762463}}},
};

template <typename F>
void expect_beta_binomial_cdfs(const char* name, int row, const F& f) {
  using stan::math::var;
  const double rel_tol = 1e-12;
  for (const auto& p : beta_binomial_cdfs_points) {
    const double* ref = p.ref[row];
    const double value_d = f(p.n, p.N, p.alpha, p.beta);
    EXPECT_NEAR(ref[0], value_d, rel_tol * std::fabs(ref[0]))
        << name << "(" << p.n << " | " << p.N << ", " << p.alpha << ", "
        << p.beta << "), double arguments";

    var alpha = p.alpha;
    var beta = p.beta;
    var value = f(p.n, p.N, alpha, beta);
    std::vector<var> vars = {alpha, beta};
    std::vector<double> grad;
    value.grad(vars, grad);
    EXPECT_NEAR(ref[0], value.val(), rel_tol * std::fabs(ref[0]))
        << name << "(" << p.n << " | " << p.N << ", " << p.alpha << ", "
        << p.beta << ")";
    EXPECT_NEAR(ref[1], grad[0], rel_tol * std::fabs(ref[1]))
        << "d/dalpha " << name << "(" << p.n << " | " << p.N << ", " << p.alpha
        << ", " << p.beta << ")";
    EXPECT_NEAR(ref[2], grad[1], rel_tol * std::fabs(ref[2]))
        << "d/dbeta " << name << "(" << p.n << " | " << p.N << ", " << p.alpha
        << ", " << p.beta << ")";
    stan::math::recover_memory();
  }
}

}  // namespace

TEST_F(AgradRev, ProbDistributionsBetaBinomial_cdf_terminating_series) {
  expect_beta_binomial_cdfs(
      "beta_binomial_cdf", 0,
      [](int n, int N, const auto& alpha, const auto& beta) {
        return stan::math::beta_binomial_cdf(n, N, alpha, beta);
      });
}

TEST_F(AgradRev, ProbDistributionsBetaBinomial_lcdf_terminating_series) {
  expect_beta_binomial_cdfs(
      "beta_binomial_lcdf", 1,
      [](int n, int N, const auto& alpha, const auto& beta) {
        return stan::math::beta_binomial_lcdf(n, N, alpha, beta);
      });
}

TEST_F(AgradRev, ProbDistributionsBetaBinomial_lccdf_terminating_series) {
  expect_beta_binomial_cdfs(
      "beta_binomial_lccdf", 2,
      [](int n, int N, const auto& alpha, const auto& beta) {
        return stan::math::beta_binomial_lccdf(n, N, alpha, beta);
      });
}

// alpha = 1 and a small cdf: beta_binomial_lcdf uses the mirrored series,
// which has a denominator parameter equal to the numerator parameter that
// ends the series
TEST_F(AgradRev, ProbDistributionsBetaBinomial_lcdf_mirror_alpha_one) {
  using stan::math::var;
  // n, N, alpha, beta, lcdf, d/dalpha, d/dbeta
  const std::vector<std::array<double, 7>> points
      = {{0, 117, 1.0, 0.1, -7.0656133635977173, -5.1910469886639266,
          9.9914602903501275},
         {0, 1000, 1.0, 0.01, -11.512935464920229, -7.4691506484692035,
          99.999000009999898}};
  for (const auto& p : points) {
    var alpha = p[2];
    var beta = p[3];
    var lcdf = stan::math::beta_binomial_lcdf(
        static_cast<int>(p[0]), static_cast<int>(p[1]), alpha, beta);
    std::vector<var> vars = {alpha, beta};
    std::vector<double> grad;
    lcdf.grad(vars, grad);
    EXPECT_NEAR(p[4], lcdf.val(), 1e-12 * std::fabs(p[4]));
    EXPECT_NEAR(p[5], grad[0], 1e-12 * std::fabs(p[5]));
    EXPECT_NEAR(p[6], grad[1], 1e-12 * std::fabs(p[6]));
    stan::math::recover_memory();
  }
}
