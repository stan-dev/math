#include <stan/math/prim.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbVonMises, boundary) {
  using stan::math::von_mises_cdf;

  EXPECT_DOUBLE_EQ(0.0, stan::math::von_mises_cdf(-stan::math::pi(), 0.6, 1.0));
  EXPECT_DOUBLE_EQ(1.0, stan::math::von_mises_cdf(stan::math::pi(), 0.6, 1.0));
  EXPECT_DOUBLE_EQ(
      1.0, stan::math::von_mises_cdf(stan::math::pi(), 6.283185307179586, 1.0));
}

TEST(ProbVonMises, throwing) {
  using stan::math::von_mises_cdf;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(stan::math::von_mises_cdf(-4, 0.6, 1.0), std::domain_error);
  EXPECT_THROW(stan::math::von_mises_cdf(4, 0.1, 1.0), std::domain_error);
  EXPECT_THROW(stan::math::von_mises_cdf(0.0, 0.1, -1.0), std::domain_error);
  EXPECT_THROW(stan::math::von_mises_cdf(nan, 0.1, 1.0), std::domain_error);
  EXPECT_THROW(stan::math::von_mises_cdf(0.0, nan, 1.0), std::domain_error);
  EXPECT_THROW(stan::math::von_mises_cdf(0.0, 0.1, nan), std::domain_error);
}

TEST(ProbVonMises, pdf_cdf_agree_test) {
  using stan::math::von_mises_cdf;
  using stan::math::von_mises_lpdf;
  double sum, t, y;
  std::vector<double> ts, ys;
  double abs_tol_int = 1e-4;

  // set random parameter values
  double mu = 0.1;
  double k = 6;

  // tabulate pdf at equispaced nodes
  int n = 1000;
  double xmin = -0.2;
  double xmax = 0.3;
  for (int i = 0; i < n; i++) {
    t = xmin + i * (xmax - xmin) / (n - 1);
    ts.push_back(t);

    y = exp(von_mises_lpdf(t, mu, k));
    ys.push_back(y);
  }

  // integrate pdf numerically with trapezoid rule
  sum = ys[0] / 2 + ys[n - 1] / 2;
  for (int i = 1; i < n - 1; i++) {
    sum += ys[i];
  }
  sum = sum / ((n - 1) / (xmax - xmin));

  // and compare to cdf
  double val2 = von_mises_cdf(xmax, mu, k) - von_mises_cdf(xmin, mu, k);

  ASSERT_NEAR(val2, sum, abs_tol_int);
}

TEST(ProbVonMises, pointwise_cdf_test) {
  using stan::math::von_mises_cdf;
  double pi = stan::math::pi();
  // check that our von_mises_cdf coincides with scipy's for mu = -pi, 0, pi
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
  ASSERT_NEAR(0.0000000001832846, von_mises_cdf_series(-1.5707463267948965, 20),
              ABS_TOL);
  ASSERT_NEAR(0.8902398103119467, von_mises_cdf_series(1.5707463267948962, 1),
              ABS_TOL);

  using stan::math::internal::von_mises_cdf_normalapprox;
  ASSERT_NEAR(0.8991261568269097,
              von_mises_cdf_normalapprox(0.16534171860998947, 60), ABS_TOL);

  using stan::math::internal::von_mises_cdf_centered;
  ASSERT_NEAR(0.1008738431730909,
              von_mises_cdf_centered(-0.16534171860998903, 60), ABS_TOL);
  ASSERT_NEAR(0.6968583086078725,
              von_mises_cdf_centered(0.16534171860998947, 10), ABS_TOL);
}
