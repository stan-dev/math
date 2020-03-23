#include <stan/math/prim.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <gtest/gtest.h>

TEST(ProbVonMises, random_cdf_test) {
  using stan::math::von_mises_cdf;
  double pi = stan::math::pi();
  // check that our von_mises_cdf coincides with scipy's
  double ABS_TOL = 1e-12;
  ASSERT_NEAR(0.0000157567385742, von_mises_cdf(-3.141492653589793, 0, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.0000157567385742,
              von_mises_cdf(-3.141492653589793, 6.283185307179586, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.0000157567385742,
              von_mises_cdf(-3.141492653589793, -6.283185307179586, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.2484164302237636, von_mises_cdf(-1.5707463267948965, 0, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.2484164302237636,
              von_mises_cdf(-1.5707463267948965, 6.283185307179586, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.2484164302237636,
              von_mises_cdf(-1.5707463267948965, -6.283185307179586, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, 0, 0.01), ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, 6.283185307179586, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, -6.283185307179586, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.7515835697762363, von_mises_cdf(1.5707463267948962, 0, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.7515835697762363,
              von_mises_cdf(1.5707463267948962, 6.283185307179586, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.7515835697762363,
              von_mises_cdf(1.5707463267948962, -6.283185307179586, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.9999842432614258, von_mises_cdf(3.141492653589793, 0, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.9999842432614258,
              von_mises_cdf(3.141492653589793, 6.283185307179586, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.9999842432614258,
              von_mises_cdf(3.141492653589793, -6.283185307179586, 0.01),
              ABS_TOL);
  ASSERT_NEAR(0.0000143649998288, von_mises_cdf(-3.141492653589793, 0, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.0000143649998288,
              von_mises_cdf(-3.141492653589793, 6.283185307179586, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.0000143649998288,
              von_mises_cdf(-3.141492653589793, -6.283185307179586, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.2341145110506426, von_mises_cdf(-1.5707463267948965, 0, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.2341145110506426,
              von_mises_cdf(-1.5707463267948965, 6.283185307179586, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.2341145110506426,
              von_mises_cdf(-1.5707463267948965, -6.283185307179586, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, 0, 0.1), ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, 6.283185307179586, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, -6.283185307179586, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.7658854889493573, von_mises_cdf(1.5707463267948962, 0, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.7658854889493573,
              von_mises_cdf(1.5707463267948962, 6.283185307179586, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.7658854889493573,
              von_mises_cdf(1.5707463267948962, -6.283185307179586, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.9999856350001711, von_mises_cdf(3.141492653589793, 0, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.9999856350001711,
              von_mises_cdf(3.141492653589793, 6.283185307179586, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.9999856350001711,
              von_mises_cdf(3.141492653589793, -6.283185307179586, 0.1),
              ABS_TOL);
  ASSERT_NEAR(0.0000046245485838, von_mises_cdf(-3.141492653589793, 0, 1),
              ABS_TOL);
  ASSERT_NEAR(0.0000046245485838,
              von_mises_cdf(-3.141492653589793, 6.283185307179586, 1), ABS_TOL);
  ASSERT_NEAR(0.0000046245485838,
              von_mises_cdf(-3.141492653589793, -6.283185307179586, 1),
              ABS_TOL);
  ASSERT_NEAR(0.1097601896880533, von_mises_cdf(-1.5707463267948965, 0, 1),
              ABS_TOL);
  ASSERT_NEAR(0.1097601896880533,
              von_mises_cdf(-1.5707463267948965, 6.283185307179586, 1),
              ABS_TOL);
  ASSERT_NEAR(0.1097601896880533,
              von_mises_cdf(-1.5707463267948965, -6.283185307179586, 1),
              ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, 0, 1), ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, 6.283185307179586, 1),
              ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, -6.283185307179586, 1),
              ABS_TOL);
  ASSERT_NEAR(0.8902398103119467, von_mises_cdf(1.5707463267948962, 0, 1),
              ABS_TOL);
  ASSERT_NEAR(0.8902398103119467,
              von_mises_cdf(1.5707463267948962, 6.283185307179586, 1), ABS_TOL);
  ASSERT_NEAR(0.8902398103119467,
              von_mises_cdf(1.5707463267948962, -6.283185307179586, 1),
              ABS_TOL);
  ASSERT_NEAR(0.9999953754514163, von_mises_cdf(3.141492653589793, 0, 1),
              ABS_TOL);
  ASSERT_NEAR(0.9999953754514163,
              von_mises_cdf(3.141492653589793, 6.283185307179586, 1), ABS_TOL);
  ASSERT_NEAR(0.9999953754514163,
              von_mises_cdf(3.141492653589793, -6.283185307179586, 1), ABS_TOL);
  ASSERT_NEAR(0.0000000000002562, von_mises_cdf(-3.141492653589793, 0, 10),
              ABS_TOL);
  ASSERT_NEAR(0.0000000000002562,
              von_mises_cdf(-3.141492653589793, 6.283185307179586, 10),
              ABS_TOL);
  ASSERT_NEAR(0.0000000000002562,
              von_mises_cdf(-3.141492653589793, -6.283185307179586, 10),
              ABS_TOL);
  ASSERT_NEAR(0.0000057188754260, von_mises_cdf(-1.5707463267948965, 0, 10),
              ABS_TOL);
  ASSERT_NEAR(0.0000057188754260,
              von_mises_cdf(-1.5707463267948965, 6.283185307179586, 10),
              ABS_TOL);
  ASSERT_NEAR(0.0000057188754260,
              von_mises_cdf(-1.5707463267948965, -6.283185307179586, 10),
              ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, 0, 10), ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, 6.283185307179586, 10),
              ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, -6.283185307179586, 10),
              ABS_TOL);
  ASSERT_NEAR(0.9999942811245741, von_mises_cdf(1.5707463267948962, 0, 10),
              ABS_TOL);
  ASSERT_NEAR(0.9999942811245741,
              von_mises_cdf(1.5707463267948962, 6.283185307179586, 10),
              ABS_TOL);
  ASSERT_NEAR(0.9999942811245741,
              von_mises_cdf(1.5707463267948962, -6.283185307179586, 10),
              ABS_TOL);
  ASSERT_NEAR(0.9999999999997439, von_mises_cdf(3.141492653589793, 0, 10),
              ABS_TOL);
  ASSERT_NEAR(0.9999999999997439,
              von_mises_cdf(3.141492653589793, 6.283185307179586, 10), ABS_TOL);
  ASSERT_NEAR(0.9999999999997439,
              von_mises_cdf(3.141492653589793, -6.283185307179586, 10),
              ABS_TOL);
  ASSERT_NEAR(0.0000000000000000, von_mises_cdf(-3.141492653589793, 0, 20),
              ABS_TOL);
  ASSERT_NEAR(0.0000000000000000,
              von_mises_cdf(-3.141492653589793, 6.283185307179586, 20),
              ABS_TOL);
  ASSERT_NEAR(0.0000000000000000,
              von_mises_cdf(-3.141492653589793, -6.283185307179586, 20),
              ABS_TOL);
  ASSERT_NEAR(0.0000000001832846, von_mises_cdf(-1.5707463267948965, 0, 20),
              ABS_TOL);
  ASSERT_NEAR(0.0000000001832846,
              von_mises_cdf(-1.5707463267948965, 6.283185307179586, 20),
              ABS_TOL);
  ASSERT_NEAR(0.0000000001832846,
              von_mises_cdf(-1.5707463267948965, -6.283185307179586, 20),
              ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, 0, 20), ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, 6.283185307179586, 20),
              ABS_TOL);
  ASSERT_NEAR(0.5000000000000000, von_mises_cdf(0.0, -6.283185307179586, 20),
              ABS_TOL);
  ASSERT_NEAR(0.9999999998167153, von_mises_cdf(1.5707463267948962, 0, 20),
              ABS_TOL);
  ASSERT_NEAR(0.9999999998167153,
              von_mises_cdf(1.5707463267948962, 6.283185307179586, 20),
              ABS_TOL);
  ASSERT_NEAR(0.9999999998167153,
              von_mises_cdf(1.5707463267948962, -6.283185307179586, 20),
              ABS_TOL);
  ASSERT_NEAR(1.0000000000000000, von_mises_cdf(3.141492653589793, 0, 20),
              ABS_TOL);
  ASSERT_NEAR(1.0000000000000000,
              von_mises_cdf(3.141492653589793, 6.283185307179586, 20), ABS_TOL);
  ASSERT_NEAR(1.0000000000000000,
              von_mises_cdf(3.141492653589793, -6.283185307179586, 20),
              ABS_TOL);
  ASSERT_NEAR(0.1008738431730909, von_mises_cdf(-0.16534171860998903, 0, 60),
	      ABS_TOL);
  ASSERT_NEAR(0.8991261568269097, von_mises_cdf(0.16534171860998947, 0, 60),
	      ABS_TOL);

  // tests for helper functions
  using stan::math::internal::von_mises_cdf_series;
  ASSERT_NEAR(0.0000000001832846,
              von_mises_cdf_series(-1.5707463267948965, 20),
	      ABS_TOL);
  ASSERT_NEAR(0.8902398103119467,
              von_mises_cdf_series(1.5707463267948962, 1),
	      ABS_TOL);

  using stan::math::internal::von_mises_cdf_normalapprox;
  ASSERT_NEAR(0.8991261568269097,
	      von_mises_cdf_normalapprox(0.16534171860998947, 60),
	      ABS_TOL);

  using stan::math::internal::von_mises_cdf_centered;
  ASSERT_NEAR(0.1008738431730909,
	      von_mises_cdf_centered(-0.16534171860998903, 60),
	      ABS_TOL);
  ASSERT_NEAR(0.6968583086078725,
	      von_mises_cdf_centered(0.16534171860998947, 10),
	      ABS_TOL);
}
