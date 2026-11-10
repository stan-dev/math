#ifndef TEST_UNIT_MATH_REV_PROB_GAMMA_LCCDF_TAIL_REFERENCE_HPP
#define TEST_UNIT_MATH_REV_PROB_GAMMA_LCCDF_TAIL_REFERENCE_HPP

#include <gtest/gtest.h>
#include <cmath>
#include <string>
#include <vector>

namespace stan {
namespace test {
namespace gamma_lccdf_tail {

/**
 * gamma_lccdf(y | alpha, beta), its first and second derivatives with
 * respect to alpha, and its derivative with respect to y.
 *
 * Computed with mpmath at 50 digits by two routes: the regularized upper
 * incomplete gamma function with numerical differentiation, and a
 * quadrature of the density scaled to 1 at the lower limit, with the
 * derivatives as ratios of positive integrals. The routes agree to 1e-47
 * where both converge; at alpha = 1e4 and 1e6 only the quadrature
 * converges.
 *
 * Rows 0 and 1 are near y = alpha + 1. Rows 2 to 9 have
 * log(Q) < log(DBL_MIN), where Q itself is subnormal or zero.
 */
struct Row {
  double alpha;
  double y;
  double beta;
  double value;
  double d_alpha;
  double d2_alpha;
  double d_y;
};

inline const std::vector<Row>& rows() {
  static const std::vector<Row> r = {
      {5, 6.5, 1, -1.4975754049608005, 5.8767262159542573e-1,
       -1.9055199558693436e-1, -4.9993873514327224e-1},
      {5, 7.5, 1, -2.0244848579180109, 6.9436601738250379e-1,
       -1.9817346506198028e-1, -5.5213820632791284e-1},
      {150, 1200, 1, -7.4345557714901795e+2, 2.0837288947176418,
       -6.6880371230448705e-3, -8.7595121899782028e-1},
      {300, 1500, 1, -7.2232706800763835e+2, 1.6119371090820549,
       -3.3382049273506081e-3, -8.0083229548002851e-1},
      {150, 600, 2, -7.4345557714901795e+2, 2.0837288947176418,
       -6.6880371230448705e-3, -1.7519024379956406},
      {10, 1000, 1, -9.5062299837015608e+2, 4.6570107370619179,
       -1.0516532158355843e-1, -9.9100906329734302e-1},
      {20, 1000, 1, -9.0807337076913739e+2, 3.9382495777101035,
       -5.1269788164528015e-2, -9.8101932786738336e-1},
      {200, 1600, 1, -9.8962693094236584e+2, 2.0826567480373627,
       -5.0120131147584039e-3, -8.757136315228946e-1},
      {1000, 3000, 1, -9.0645450699007447e+2, 1.0996116242418308,
       -1.0002512857967094e-3, -6.6716616809812125e-1},
      {1e6, 1.1e6, 1, -4.6953444144485314e+3, 9.5320677605584263e-2,
       -9.9990056594105228e-7, -9.0919088810213022e-2},
      {700, 2100, 1, -6.3586021219348883e+2, 1.1000395060722886,
       -1.4290853740612502e-3, -6.6737993613863034e-1},
      {1e4, 10100, 1, -1.84104688482372, 1.5226403888940706e-2,
       -8.0334369099635953e-5, -1.5150648133166995e-2},
      {200, 200, 1, -7.1213311417597682e-1, 5.7524325573397242e-2,
       -3.3091679741037347e-3, -5.7476404689539276e-2},
  };
  return r;
}

/** Tolerance for second derivatives: the continued fraction stops on the
 * convergence of the value and the first derivative only. */
inline double d2_tol(const Row& r) { return r.y < r.alpha + 3 ? 5e-8 : 1e-9; }

inline std::string label(const Row& r) {
  return "alpha=" + std::to_string(r.alpha) + " y=" + std::to_string(r.y)
         + " beta=" + std::to_string(r.beta);
}

inline void expect_rel(double expected, double actual, double tol,
                       const std::string& what) {
  ASSERT_TRUE(std::isfinite(actual)) << what << ": not finite, " << actual;
  const double err = std::fabs(actual - expected) / std::fabs(expected);
  EXPECT_LT(err, tol) << what << ": relative error " << err << ", expected "
                      << expected << ", got " << actual;
}

}  // namespace gamma_lccdf_tail
}  // namespace test
}  // namespace stan

#endif
