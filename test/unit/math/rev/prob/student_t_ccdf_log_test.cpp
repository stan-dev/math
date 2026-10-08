#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

namespace {

// student_t_lccdf and its nu and sigma gradients against a high-precision
// reference (mpmath, 50 digits, the tail as 0.5 I_{nu/(nu+t^2)}(nu/2, 1/2)
// differentiated numerically).
void expect_student_t_lccdf(double y, double nu_in, double sigma_in,
                            double expected, double dnu, double dsigma,
                            double rtol) {
  using stan::math::var;
  var nu = nu_in;
  var sigma = sigma_in;
  var lp = stan::math::student_t_lccdf(y, nu, 0.0, sigma);
  std::vector<var> x{nu, sigma};
  std::vector<double> g;
  lp.grad(x, g);
  EXPECT_TRUE(std::isfinite(lp.val())) << "y=" << y << " nu=" << nu_in;
  EXPECT_TRUE(std::isfinite(g[0]) && std::isfinite(g[1]))
      << "y=" << y << " nu=" << nu_in;
  EXPECT_NEAR(expected, lp.val(), rtol * std::fabs(expected));
  EXPECT_NEAR(dnu, g[0], rtol * std::fabs(dnu));
  EXPECT_NEAR(dsigma, g[1], rtol * std::fabs(dsigma));
}

// f(y, nu, sigma) with location 0: the value, then the y, nu and sigma
// gradients, against references (mpmath, 60 digits, two routes)
template <typename F>
void expect_student_t_fn(F f, double y_in, double nu_in, double sigma_in,
                         const std::vector<double>& expected, double rtol) {
  using stan::math::var;
  var y = y_in;
  var nu = nu_in;
  var sigma = sigma_in;
  var p = f(y, nu, sigma);
  std::vector<var> x{y, nu, sigma};
  std::vector<double> g;
  p.grad(x, g);
  EXPECT_NEAR(expected[0], p.val(), rtol * std::fabs(expected[0]))
      << "y=" << y_in << " nu=" << nu_in;
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(expected[i + 1], g[i], rtol * std::fabs(expected[i + 1]))
        << "gradient " << i << " y=" << y_in << " nu=" << nu_in;
  }
}

const auto student_t_cdf_fn = [](auto& y, auto& nu, auto& sigma) {
  return stan::math::student_t_cdf(y, nu, 0.0, sigma);
};
const auto student_t_lcdf_fn = [](auto& y, auto& nu, auto& sigma) {
  return stan::math::student_t_lcdf(y, nu, 0.0, sigma);
};
const auto student_t_lccdf_fn = [](auto& y, auto& nu, auto& sigma) {
  return stan::math::student_t_lccdf(y, nu, 0.0, sigma);
};

}  // namespace

TEST_F(AgradRev, ProbDistributionsStudentT_lccdf_deep_tail) {
  // q = nu / t^2 >= 2, the branch that formed the complement in linear
  // space; the true value is -50.8, below eps = e^-36.7
  expect_student_t_lccdf(10, 1000, 1, -50.838951454, -0.0022457674682,
                         91.799148523, 1e-9);
  expect_student_t_lccdf(8, 200, 1, -30.638906854, -0.018204242375,
                         49.213616295, 1e-9);
  // q < 2, the other branch, unchanged
  expect_student_t_lccdf(20, 50, 1, -57.754062809, -0.66295856977, 44.550802193,
                         1e-9);
  expect_student_t_lccdf(50, 1000, 1, -630.58671311, -0.26959552076,
                         714.57063175, 1e-9);
  expect_student_t_lccdf(4, 4, 1, -4.8202159939, -0.49244669515, 3.3270509831,
                         1e-9);
}

TEST_F(AgradRev, ProbDistributionsStudentT_lccdf_large_t) {
  // stan-dev/math#2923: q = nu / t^2 < 2, where 1 - r with r = 1 / (1 + q)
  // was formed by subtraction. It loses digits from |t| = 10 on and is 0
  // once q is below eps, which made the tail probability 0.
  expect_student_t_lccdf(1e8, 0.1, 1, -2.7158255239161981, -19.691030420120541,
                         0.10000000000000001, 1e-12);
  expect_student_t_lccdf(1e29, 0.1, 1, -7.5512542192036944, -68.045317372995498,
                         0.10000000000000001, 1e-12);
  // Cauchy, t = 1e10 with sigma = 2.5
  expect_student_t_lccdf(2.5e10, 1, 2.5, -24.170580815789858,
                         -22.83270374938051, 0.40000000000000002, 1e-12);
  // finite before the change, but with a relative error of 3e-5
  expect_student_t_lccdf(1e7, 4, 1, -63.373770315165238, -15.034762317625022,
                         3.9999999999998668, 1e-12);
}

TEST_F(AgradRev, ProbDistributionsStudentT_cdf_lcdf_large_t) {
  // stan-dev/math#2923: the reported case nu = 0.1 at y = -1e8, where the
  // cdf was exactly 0, its mirror at y = 1e8, where it was exactly 1, and
  // the Cauchy lcdf at y = -1e10. Value, then the y, nu and sigma gradients.
  expect_student_t_fn(student_t_cdf_fn, -1e8, 0.1, 1.0,
                      {0.066150321787786709, 6.6150321787786714e-11,
                       -1.3025679986240706, 0.0066150321787786714},
                      1e-12);
  expect_student_t_fn(student_t_cdf_fn, 1e8, 0.1, 1.0,
                      {0.93384967821221332, 6.6150321787786714e-11,
                       1.3025679986240706, -0.0066150321787786714},
                      1e-12);
  expect_student_t_fn(student_t_lcdf_fn, -1e10, 1.0, 1.0,
                      {-24.170580815789858, 1e-10, -22.83270374938051, 1.0},
                      1e-12);
}

TEST_F(AgradRev, ProbDistributionsStudentT_cdf_lccdf_near_zero) {
  // q >= 2 with |t| small. The complement 1 - I_r(1/2, nu/2) and its nu
  // gradient must be evaluated from r, not from 1 - r, which is close to 1
  // here: from 1 - r the nu gradient had a relative error of 1e-3 at
  // t = 1e-7. nu = 1 is Boost's arcsine case a = b = 1/2.
  expect_student_t_fn(student_t_cdf_fn, -1e-7, 1.0, 1.0,
                      {0.49999996816901138, 0.31830988618378747,
                       -6.1480657060756583e-09, 3.1830988618378744e-08},
                      1e-11);
  expect_student_t_fn(student_t_lccdf_fn, 1e-3, 1.0, 1.0,
                      {-0.69378400284838371, -0.63702467811884134,
                       -0.00012303970872281812, 0.00063702467811884139},
                      1e-11);
  expect_student_t_fn(student_t_cdf_fn, -1e-5, 30.0, 1.0,
                      {0.49999604367815115, 0.39563218487365676,
                       -1.0983690983393141e-09, 3.9563218487365682e-06},
                      1e-11);
}

TEST_F(AgradRev, ProbDistributionsStudentT_lccdf_branch_seam) {
  // the two branches meet at q = 2; the value and the gradients must be
  // continuous across it
  using stan::math::var;
  double v[3], d[3];
  int k = 0;
  for (double y : {1.9999, 2.0, 2.0001}) {
    var nu = 8;
    var lp = stan::math::student_t_lccdf(y, nu, 0.0, 1.0);
    lp.grad();
    v[k] = lp.val();
    d[k] = nu.adj();
    ++k;
  }
  EXPECT_NEAR(v[1], 0.5 * (v[0] + v[2]), 1e-8);
  EXPECT_NEAR(d[1], 0.5 * (d[0] + d[2]), 1e-8);
}

TEST_F(AgradRev, ProbDistributionsStudentT_lcdf_cdf_deep_tail) {
  using stan::math::var;
  {
    var nu = 1000;
    var sigma = 1;
    var lp = stan::math::student_t_lcdf(-10.0, nu, 0.0, sigma);
    lp.grad();
    EXPECT_TRUE(std::isfinite(lp.val()));
    EXPECT_NEAR(-50.838951454, lp.val(), 1e-9 * 50.9);
    EXPECT_NEAR(-0.0022457674682, nu.adj(), 1e-9 * 0.00225);
  }
  {
    // the cdf itself: the value is representable and the gradients were NaN
    var nu = 1000;
    var sigma = 1;
    var p = stan::math::student_t_cdf(-10.0, nu, 0.0, sigma);
    p.grad();
    EXPECT_TRUE(std::isfinite(nu.adj()) && std::isfinite(sigma.adj()));
    EXPECT_NEAR(8.33535e-23, p.val(), 1e-5 * 8.3e-23);
  }
}
