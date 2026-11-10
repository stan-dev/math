#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/rev/prob/poisson_log_cdfs_reference.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <string>

namespace {

/**
 * Relative error; a reference below 1e-300 (zero or subnormal, with few
 * significant bits) needs an absolute error below 1e-300.
 */
void expect_close(double expected, double actual, double tol,
                  const std::string& what) {
  ASSERT_TRUE(std::isfinite(actual)) << what << ": not finite, " << actual;
  if (std::fabs(expected) < 1e-300) {
    EXPECT_LE(std::fabs(actual - expected), 1e-300)
        << what << ": expected " << expected << ", got " << actual;
    return;
  }
  const double err = std::fabs(actual - expected) / std::fabs(expected);
  EXPECT_LT(err, tol) << what << ": relative error " << err << ", expected "
                      << expected << ", got " << actual;
}

}  // namespace

TEST_F(AgradRev, ProbDistributionsPoissonLogCdfs_reference) {
  using stan::math::var;
  namespace ref = stan::test::poisson_log_cdfs;
  for (const auto& r : ref::rows()) {
    const int n = static_cast<int>(r.y);
    const std::string at = r.function + " at n = " + std::to_string(n)
                           + ", alpha = " + std::to_string(r.p1);
    ASSERT_EQ(1, r.n_grad) << at;
    var alpha = r.p1;
    var lp;
    if (r.function == "poisson_log_cdf") {
      lp = stan::math::poisson_log_cdf(n, alpha);
    } else if (r.function == "poisson_log_lcdf") {
      lp = stan::math::poisson_log_lcdf(n, alpha);
    } else {
      ASSERT_EQ("poisson_log_lccdf", r.function);
      lp = stan::math::poisson_log_lccdf(n, alpha);
    }
    lp.grad();
    // exp(alpha) is rounded, and near the mean of a large count the value
    // changes by thousands of times the relative change of the rate
    expect_close(r.value, lp.val(), 1e-11, "value of " + at);
    expect_close(r.grad[0], alpha.adj(), 1e-9, "gradient of " + at);
    stan::math::recover_memory();
  }
}
