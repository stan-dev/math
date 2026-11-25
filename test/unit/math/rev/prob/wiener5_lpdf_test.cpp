#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>

#include <gtest/gtest.h>
#include <cmath>
#include <vector>

TEST_F(AgradRev, ProbWiener5_small_density) {
  /* Points where the density is far below 1. The series tolerance of the
   * density and of the gradients must be relative to the density there;
   * with an absolute tolerance the series stop after one or two terms and
   * the log density is wrong by up to 0.47 (point 2).
   * References: the large-time series of the first-passage time density
   * (400 terms, 50 digits, mpmath), with the drift variability in closed
   * form (Blurton, Kesselmeier and Gondan, 2017):
   *   p(t) = exp((sv^2 a^2 w^2 - 2 a v w - v^2 t) / (2 (1 + sv^2 t)))
   *          / sqrt(1 + sv^2 t) * (pi / a^2)
   *          * sum_k k sin(k pi w) exp(-k^2 pi^2 t / (2 a^2))
   * for the lower boundary; the upper boundary uses 1 - w and -v.
   * Gradients by mpmath.diff.
   */
  using stan::math::var;
  // y, a, w, v, sv (t0 = 0)
  const std::vector<std::vector<double>> points
      = {{1.0, 2.0, 0.2, -6.0, 0.0},
         {10.0, 8.0, 1e-4, -10.0, 0.0},
         {3.0, 3.0, 0.1, -4.0, 0.0},
         {1.0, 2.0, 0.2, -6.0, 0.5},
         {1.40338, 2.99164, 0.217078, -4.91775, 0.0470162}};
  // value, then d/dy, d/da, d/dt0, d/dw, d/dv, d/dsv
  const std::vector<std::vector<double>> refs
      = {{-29.689679323384301, -18.914999540875172, -4.8850004591248286,
          18.914999540875172, 16.520222792214895, 7.5999999999999996, 0},
         {-592.31616885779931, -49.948518526705826, -10.377703683235437,
          49.948518526705826, 10080.000219654299, 107.9992, 0},
         {-38.699512995035121, -8.5021097351560204, -3.2624471963546275,
          8.5021097351560204, 21.697496903138553, 14.699999999999999, 0},
         {-24.025251099041405, -11.050199540875171, -3.6690004591248289,
          11.050199540875171, 13.480222792214894, 6.0800000000000001,
          18.083200000000001},
         {-31.029062778875129, -11.882570552710849, -4.6066508752548998,
          11.882570552710849, 19.518498453098601, 9.2151056308643184,
          3.9267522089987508}};
  for (size_t i = 0; i < points.size(); ++i) {
    const auto& p = points[i];
    std::vector<var> x{p[0], p[1], 0.0, p[2], p[3], p[4]};
    var lp = stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4], x[5]);
    lp.grad();
    EXPECT_NEAR(lp.val(), refs[i][0], 1e-6) << "point " << i;
    for (size_t k = 0; k < x.size(); ++k) {
      EXPECT_NEAR(x[k].adj(), refs[i][k + 1],
                  1e-8 * std::fmax(1.0, std::fabs(refs[i][k + 1])))
          << "point " << i << ", argument " << k;
    }
    stan::math::recover_memory();
  }
}
