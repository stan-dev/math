#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <vector>

TEST_F(AgradRev, ProbDistributionsBeta_lcdf_large_beta_d_y) {
  using stan::math::var;
  // Reference values from mpmath at 50 digits. Columns: y, alpha, beta,
  // beta_lcdf(y | alpha, beta), its derivative with respect to y. y and
  // 1 - y are exact, so the derivative depends on beta(alpha, beta).
  const std::vector<std::array<double, 5>> rows{{
      {1.862645149230957e-09, 2.0, 1073741824.0, -0.5208858063911924,
       489281430.6449945},
      {2.3283064365386963e-10, 2.0, 2147483648.0, -2.40568139077329,
       7219826854.641048},
  }};
  // The lgamma-difference formula for beta() gives relative derivative
  // errors of 2e-6 to 3e-6 at these points.  Windows x64 (MinGW) gives
  // relative errors up to about 2e-8 in the value and the derivative.
  const double tol = 1e-7;
  for (const auto& r : rows) {
    SCOPED_TRACE(::testing::Message() << "beta = " << r[2]);
    var y = r[0];
    var lp = stan::math::beta_lcdf(y, r[1], r[2]);
    lp.grad();
    EXPECT_NEAR(r[3], lp.val(), tol * std::fabs(r[3]));
    EXPECT_NEAR(r[4], y.adj(), tol * r[4]);
    stan::math::recover_memory();
  }
}
