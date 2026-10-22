#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>

#include <gtest/gtest.h>
#include <cmath>
#include <vector>

TEST_F(AgradRev, ProbWiener4Lcdf_small_time_gradients) {
  /* Gradients of the 4-parameter lcdf at points where the small-time series
   * is used. Its derivatives with respect to a and w had sign errors (errors
   * of 1.4e-4 in d/da and 3.8e-4 in d/dw at the first point). References:
   * mpmath (50 digits): F = P - Q, P in closed form and Q = P - F by its
   * large-time series (199 terms); gradients by mpmath.diff.
   */
  using stan::math::var;
  // y, a, t0, w, v
  const std::vector<std::vector<double>> points
      = {{0.4693881426215285, 2.2901173915641948, 0.1, 0.35423574118263701,
          4.2329665942670687},
         {0.25358060660214943, 1.8799221079023491, 0.1, 0.3028584205531582,
          -0.69503918555849786}};
  // value, then d/dy, d/da, d/dt0, d/dw, d/dv
  const std::vector<std::vector<double>> refs
      = {{-0.45962394210453974, 4.1206910657932813, -0.67059830839645508,
          -4.1206910657932813, 2.3789919951513356, 0.35637417773431785},
         {-8.043108690545884, 39.059443949324603, -6.9008508266698447,
          -39.059443949324603, 18.608978302067688, 1.403540781731772}};
  for (size_t i = 0; i < points.size(); ++i) {
    std::vector<var> x(points[i].begin(), points[i].end());
    var lp = stan::math::wiener_lcdf_unnorm(x[0], x[1], x[2], x[3], x[4]);
    lp.grad();
    EXPECT_NEAR(lp.val(), refs[i][0], 1e-10) << "point " << i;
    for (size_t k = 0; k < x.size(); ++k) {
      EXPECT_NEAR(x[k].adj(), refs[i][k + 1],
                  1e-8 * std::fmax(1.0, std::fabs(refs[i][k + 1])))
          << "point " << i << ", argument " << k;
    }
    stan::math::recover_memory();
  }
}
