#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>

#include <gtest/gtest.h>
#include <cmath>
#include <vector>

TEST_F(AgradRev, ProbWiener4Lccdf_complement) {
  /* The ccdf is P - F, P the probability of the upper boundary. Where F is
   * close to P, the difference P - F loses digits, down to -inf for a log
   * ccdf below about -38; the large-time series of P - F has no
   * cancellation. F / P at the points: 1, 1, 0.82, 0.83, 0.13, 0.079.
   * References: P in closed form and the large-time series of P - F
   * (199 terms, 50 digits, mpmath); gradients by mpmath.diff.
   */
  using stan::math::var;
  // y, a, t0, w, v
  const std::vector<std::vector<double>> points
      = {{2.1, 0.55, 0.1, 0.5, -3.5}, {3.0, 1.0, 0.1, 0.6, 4.5},
         {1.5, 2.0, 0.1, 0.3, 1.0},   {0.5, 1.0, 0.1, 0.5, 0.5},
         {0.25, 1.5, 0.1, 0.5, 1.0},  {0.4, 2.0, 0.1, 0.4, -0.5}};
  // value, then d/dy, d/da, d/dt0, d/dw, d/dv
  const std::vector<std::vector<double>> refs
      = {{-46.609661106618847, -22.438395704279932, 115.9002611058362,
          22.438395704279932, 1.925, 7.4309826311170903},
         {-43.490907373843304, -15.059802200544679, 29.077213588597631,
          15.059802200544679, -5.5207653306919257, -12.948808705458113},
         {-2.0351620810484594, -1.7180376564123054, 2.1160210729321483,
          1.7180376564123054, 0.30079638170035061, -0.57846329209015851},
         {-2.2505184333541364, -5.0598003351182337, 4.1484312199614369,
          5.0598003351182337, -0.49142007849811765, 0.20118190373369988},
         {-0.33998143949291396, -2.179378824828539, 0.58202009464491435,
          2.179378824828539, -0.12772855946593473, 0.21921649451880984},
         {-1.7336108762577432, -0.7908861264789042, -0.20155913011329556,
          0.7908861264789042, 3.0587802685897929, 1.7552998722278673}};
  // Relative gradient tolerance: the series (F >= P / 2) is accurate to
  // rounding; where F < P / 2 the gradients of F enter, which follow the
  // default precision_derivatives = 1e-4.
  const std::vector<double> grad_tol = {1e-8, 1e-8, 1e-8, 1e-8, 1e-4, 1e-4};
  for (size_t i = 0; i < points.size(); ++i) {
    std::vector<var> x(points[i].begin(), points[i].end());
    var lp = stan::math::wiener_lccdf_unnorm(x[0], x[1], x[2], x[3], x[4]);
    lp.grad();
    EXPECT_NEAR(lp.val(), refs[i][0], 1e-9) << "point " << i;
    for (size_t k = 0; k < x.size(); ++k) {
      EXPECT_NEAR(x[k].adj(), refs[i][k + 1],
                  grad_tol[i] * std::fmax(1.0, std::fabs(refs[i][k + 1])))
          << "point " << i << ", argument " << k;
    }
    stan::math::recover_memory();
  }
}
