#ifndef TEST_UNIT_MATH_REV_PROB_SKEW_NORMAL_CDF_TEST_HELPERS_HPP
#define TEST_UNIT_MATH_REV_PROB_SKEW_NORMAL_CDF_TEST_HELPERS_HPP

#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

namespace skew_normal_cdf_test {

// One reference point: the arguments, the log cdf or log ccdf, and its
// partials with respect to y, mu, sigma and alpha.
struct ref_row {
  double y, mu, sigma, alpha;
  double value, d_y, d_mu, d_sigma, d_alpha;
};

// Compares the value and all four partials of f(y, mu, sigma, alpha) with
// the reference, to the relative tolerance rel_tol.
template <typename F>
void expect_row(const ref_row& r, const F& f, double rel_tol) {
  using stan::math::var;
  var y = r.y;
  var mu = r.mu;
  var sigma = r.sigma;
  var alpha = r.alpha;
  var lp = f(y, mu, sigma, alpha);
  std::vector<var> x{y, mu, sigma, alpha};
  std::vector<double> g;
  lp.grad(x, g);
  const double ref[5] = {r.value, r.d_y, r.d_mu, r.d_sigma, r.d_alpha};
  const double got[5] = {lp.val(), g[0], g[1], g[2], g[3]};
  const char* name[5] = {"value", "d_y", "d_mu", "d_sigma", "d_alpha"};
  for (int i = 0; i < 5; ++i) {
    EXPECT_NEAR(got[i], ref[i], rel_tol * std::fabs(ref[i]))
        << name[i] << " at y = " << r.y << ", mu = " << r.mu
        << ", sigma = " << r.sigma << ", alpha = " << r.alpha;
  }
  stan::math::recover_memory();
}

}  // namespace skew_normal_cdf_test

#endif
