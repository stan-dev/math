#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/rev/prob/skew_normal_cdf_test_helpers.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace skew_normal_lccdf_test {

using skew_normal_cdf_test::ref_row;

// References: mpmath at 40 digits. The tail probability is computed twice,
// from the Owen's T integral without subtraction and from the integral of
// the density over the tail; the two agree to at least 38 digits.
// A partial below the double range is written as 0.

// Upper tail. Computing the tail as 1 - Phi(z) loses all digits from
// z = 8 on and gives -inf from z = 9 on when alpha <= 0.
// clang-format off
// NOLINTBEGIN(whitespace/line_length)
const std::vector<ref_row> upper_tail = {
  {9.0, 0.0, 1.0, 0.0, -43.628149113332115, -9.1085231050028688, 9.1085231050028688, 81.976707945025819, 7.2675499571979654},
  {20.0, 0.0, 1.0, 0.0, -203.91715537109726, -20.049753068527851, 20.049753068527851, 400.99506137055701, 15.997388421288246},
  {8.0, 0.0, 1.0, 1.0, -34.320289979354605, -8.1213681122361102, 8.1213681122361102, 64.970944897888881, 2.0515676636104562e-14},
  {30.0, 0.0, 1.0, 5.0, -453.62809677578325, -30.033259667433677, 30.033259667433677, 900.99779002301031, 0.0},
  {12.0, 0.0, 1.0, 1000.0, -74.717525821008851, -12.082214175254284, 12.082214175254284, 144.98657010305141, 0.0},
  {10.0, 0.0, 1.0, -0.1, -54.394085623496543, -10.249947932778125, 10.249947932778125, 102.49947932778125, 15.477779377320044},
  {25.0, 3.0, 2.0, 2.0, -63.13178691386377, -5.5447325148575862, 5.5447325148575862, 60.992057663433448, 7.040216854905144e-106},
  {35.0, 0.0, 1.0, 0.1, -616.28218427221583, -35.028438615296751, 35.028438615296751, 1225.9953515353863, 0.030273094099071228},
};

// Short tail of the ccdf, z > 0 and alpha < 0: 1 - Phi(z) + 2 T(z, alpha)
// cancels there, for alpha |z| >= 6 to all digits.
const std::vector<ref_row> short_tail = {
  {1.0, 0.0, 1.0, -10.0, -58.591263384976217, -102.94417817175902, 102.94417817175902, 102.94417817175902, 10.292474347248229},
  {1.5, 0.0, 1.0, -5.0, -36.122320271091644, -40.271258082380923, 40.271258082380923, 60.406887123571384, 11.816464401459614},
  {40.0, 0.0, 1.0, -0.1, -814.28589778635187, -40.44750446065164, 40.44750446065164, 1617.9001784260656, 169.22303349079141},
  {6.0, 1.0, 0.5, -3.0, -509.15419733250449, -200.39876451091412, 200.39876451091412, 2003.9876451091412, 300.93140657945882},
};

// Bulk and the side where the ccdf is close to 1.
const std::vector<ref_row> bulk = {
  {0.5, 0.0, 1.0, -3.0, -5.0562417512463643, -7.3854067902549387, 7.3854067902549387, 3.6927033951274694, 1.4317919510486336},
  {-2.0, 0.0, 1.0, 1.0, -5.1770248847042769e-4, -0.0024578753432330854, 0.0024578753432330854, -0.0049157506864661707, 0.0029165339711515148},
  {1.3, -0.4, 0.7, 2.5, -4.189197880587666, -3.9396893672846619, 3.9396893672846619, 9.5678170348341799, 1.5017488395285856e-9},
};
// NOLINTEND
// clang-format on

constexpr double rel_tol = 1e-12;

struct lccdf_functor {
  template <typename T1, typename T2, typename T3, typename T4>
  auto operator()(const T1& y, const T2& mu, const T3& sigma,
                  const T4& alpha) const {
    return stan::math::skew_normal_lccdf(y, mu, sigma, alpha);
  }
};

}  // namespace skew_normal_lccdf_test

TEST_F(AgradRev, ProbDistributionsSkewNormal_lccdf_upper_tail) {
  using skew_normal_lccdf_test::lccdf_functor;
  for (const auto& r : skew_normal_lccdf_test::upper_tail) {
    skew_normal_cdf_test::expect_row(r, lccdf_functor(),
                                     skew_normal_lccdf_test::rel_tol);
  }
}

TEST_F(AgradRev, ProbDistributionsSkewNormal_lccdf_short_tail) {
  using skew_normal_lccdf_test::lccdf_functor;
  for (const auto& r : skew_normal_lccdf_test::short_tail) {
    skew_normal_cdf_test::expect_row(r, lccdf_functor(),
                                     skew_normal_lccdf_test::rel_tol);
  }
}

TEST_F(AgradRev, ProbDistributionsSkewNormal_lccdf_bulk) {
  using skew_normal_lccdf_test::lccdf_functor;
  for (const auto& r : skew_normal_lccdf_test::bulk) {
    skew_normal_cdf_test::expect_row(r, lccdf_functor(),
                                     skew_normal_lccdf_test::rel_tol);
  }
}

TEST_F(AgradRev, ProbDistributionsSkewNormal_lccdf_vectorized_tail) {
  using skew_normal_lccdf_test::rel_tol;
  using skew_normal_lccdf_test::upper_tail;
  using stan::math::var;
  // the first three rows have mu = 0 and sigma = 1
  std::vector<var> y{upper_tail[0].y, upper_tail[1].y, upper_tail[2].y};
  std::vector<var> alpha{upper_tail[0].alpha, upper_tail[1].alpha,
                         upper_tail[2].alpha};
  var mu = 0.0;
  var sigma = 1.0;
  var lp = stan::math::skew_normal_lccdf(y, mu, sigma, alpha);
  std::vector<var> x{y[0], y[1], y[2], alpha[0], alpha[1], alpha[2], mu, sigma};
  std::vector<double> g;
  lp.grad(x, g);
  double value = 0;
  double d_mu = 0;
  double d_sigma = 0;
  for (int i = 0; i < 3; ++i) {
    value += upper_tail[i].value;
    d_mu += upper_tail[i].d_mu;
    d_sigma += upper_tail[i].d_sigma;
    EXPECT_NEAR(g[i], upper_tail[i].d_y,
                rel_tol * std::fabs(upper_tail[i].d_y));
    EXPECT_NEAR(g[3 + i], upper_tail[i].d_alpha,
                rel_tol * std::fabs(upper_tail[i].d_alpha));
  }
  EXPECT_NEAR(lp.val(), value, rel_tol * std::fabs(value));
  EXPECT_NEAR(g[6], d_mu, rel_tol * std::fabs(d_mu));
  EXPECT_NEAR(g[7], d_sigma, rel_tol * std::fabs(d_sigma));
}

TEST_F(AgradRev, ProbDistributionsSkewNormal_lccdf_error_names_function) {
  EXPECT_THROW_MSG(stan::math::skew_normal_lccdf(1.0, 0.0, -1.0, 0.0),
                   std::domain_error, "skew_normal_lccdf");
  EXPECT_THROW_MSG(stan::math::skew_normal_lcdf(1.0, 0.0, -1.0, 0.0),
                   std::domain_error, "skew_normal_lcdf");
}
