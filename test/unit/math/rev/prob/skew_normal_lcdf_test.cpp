#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/rev/prob/skew_normal_cdf_test_helpers.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

namespace skew_normal_lcdf_test {

using skew_normal_cdf_test::ref_row;

// References: mpmath at 40 digits, computed from two integral
// representations of the tail probability without subtraction; the two
// agree to at least 38 digits.

// Short tail, z < 0 and alpha > 0: Phi(z) - 2 T(z, alpha) cancels there,
// for alpha |z| >= 6 to all digits, and Phi(z) underflows below z = -38.5.
// clang-format off
// NOLINTBEGIN(whitespace/line_length)
const std::vector<ref_row> short_tail = {
  {-1.0, 0.0, 1.0, 10.0, -58.591263384976217, 102.94417817175902, -102.94417817175902, 102.94417817175902, -10.292474347248229},
  {-1.5, 0.0, 1.0, 5.0, -36.122320271091644, 40.271258082380923, -40.271258082380923, 60.406887123571384, -11.816464401459614},
  {-4.0, 0.0, 1.0, 2.0, -46.258745719683616, 20.48140102231737, -20.48140102231737, 81.925604089269479, -33.267399431313681},
  {-0.3, 0.0, 1.0, 20.0, -25.847782710226922, 126.49748260419822, -126.49748260419822, 37.949244781259463, -1.9427245539567025},
  {-8.0, 0.0, 1.0, 1.0, -70.0268743198291, 16.242736224472225, -16.242736224472225, 129.9418897957778, -65.956620014445561},
  {-30.0, 0.0, 1.0, 0.5, -569.98328250194864, 37.566259867497551, -37.566259867497551, 1126.9877960249265, -452.78122634853445},
  {-40.0, 0.0, 1.0, 0.1, -814.28589778635187, 40.44750446065164, -40.44750446065164, 1617.9001784260656, -169.22303349079141},
  {-0.001, 0.0, 1.0, 10000.0, -64.989254270556036, 1.0194383133447667e+5, -1.0194383133447667e+5, 101.94383133447667, -0.010294383031485114},
  {-1.0, 2.0, 3.0, 10.0, -58.591263384976217, 34.314726057253007, -34.314726057253007, 34.314726057253007, -10.292474347248229},
};

// alpha > 1 outside the short tail. In the first row, close to 0 with a
// large alpha, Phi(z) - 2 T(z, alpha) loses about 3 digits.
const std::vector<ref_row> large_alpha = {
  {1e-08, 0.0, 1.0, 1000.0, -8.0524729650477151, 1253.3088471252344, -1253.3088471252344, -1.2533088471252344e-5, -9.9998680025379521e-4},
  {-0.5, 0.0, 1.0, 4.0, -6.5537546798400579, 11.243374181650123, -11.243374181650123, 5.6216870908250613, -1.5695853087783158},
  {0.7, 0.0, 1.0, 2.5, -0.65454488357280261, 1.153579725920523, -1.153579725920523, -0.80750580814436608, -0.014300847594313707},
};

// alpha <= 1 outside the short tail.
const std::vector<ref_row> small_alpha = {
  {-2.0, 0.0, 1.0, 0.5, -5.2030115545466701, 3.1149046969528357, -3.1149046969528357, 6.2298093939056714, -3.8005208281618404},
  {-0.5, 0.0, 1.0, -2.0, -0.52475764856158949, 1.0012132522892172, -1.0012132522892172, 0.50060662614460858, -0.057589780452459821},
  {1.2, 0.3, 0.8, 0.9, -0.272577298571282, 0.58738745536999308, -0.58738745536999308, -0.66081088729124216, -0.073470592302399931},
};
// NOLINTEND
// clang-format on

constexpr double rel_tol = 1e-12;

struct lcdf_functor {
  template <typename T1, typename T2, typename T3, typename T4>
  auto operator()(const T1& y, const T2& mu, const T3& sigma,
                  const T4& alpha) const {
    return stan::math::skew_normal_lcdf(y, mu, sigma, alpha);
  }
};

void expect_rows(const std::vector<ref_row>& rows) {
  for (const auto& r : rows) {
    skew_normal_cdf_test::expect_row(r, lcdf_functor(), rel_tol);
  }
}

}  // namespace skew_normal_lcdf_test

TEST_F(AgradRev, ProbDistributionsSkewNormal_lcdf_short_tail) {
  skew_normal_lcdf_test::expect_rows(skew_normal_lcdf_test::short_tail);
}

TEST_F(AgradRev, ProbDistributionsSkewNormal_lcdf_large_alpha) {
  skew_normal_lcdf_test::expect_rows(skew_normal_lcdf_test::large_alpha);
}

TEST_F(AgradRev, ProbDistributionsSkewNormal_lcdf_small_alpha) {
  skew_normal_lcdf_test::expect_rows(skew_normal_lcdf_test::small_alpha);
}

TEST_F(AgradRev, ProbDistributionsSkewNormal_lcdf_vectorized_branches) {
  using skew_normal_lcdf_test::rel_tol;
  using stan::math::var;
  // one row from each representation, all with mu = 0 and sigma = 1
  const skew_normal_cdf_test::ref_row rows[3]
      = {skew_normal_lcdf_test::short_tail[0],
         skew_normal_lcdf_test::large_alpha[1],
         skew_normal_lcdf_test::small_alpha[0]};
  std::vector<var> y{rows[0].y, rows[1].y, rows[2].y};
  std::vector<var> alpha{rows[0].alpha, rows[1].alpha, rows[2].alpha};
  var mu = 0.0;
  var sigma = 1.0;
  var lp = stan::math::skew_normal_lcdf(y, mu, sigma, alpha);
  std::vector<var> x{y[0], y[1], y[2], alpha[0], alpha[1], alpha[2], mu, sigma};
  std::vector<double> g;
  lp.grad(x, g);
  double value = 0;
  double d_mu = 0;
  double d_sigma = 0;
  for (int i = 0; i < 3; ++i) {
    value += rows[i].value;
    d_mu += rows[i].d_mu;
    d_sigma += rows[i].d_sigma;
    EXPECT_NEAR(g[i], rows[i].d_y, rel_tol * std::fabs(rows[i].d_y));
    EXPECT_NEAR(g[3 + i], rows[i].d_alpha,
                rel_tol * std::fabs(rows[i].d_alpha));
  }
  EXPECT_NEAR(lp.val(), value, rel_tol * std::fabs(value));
  EXPECT_NEAR(g[6], d_mu, rel_tol * std::fabs(d_mu));
  EXPECT_NEAR(g[7], d_sigma, rel_tol * std::fabs(d_sigma));
}

TEST_F(AgradRev, ProbDistributionsSkewNormal_lcdf_infinite_y) {
  using stan::math::var;
  var y = -std::numeric_limits<double>::infinity();
  var alpha = 2.0;
  var lp = stan::math::skew_normal_lcdf(y, 0.0, 1.0, alpha);
  EXPECT_EQ(lp.val(), -std::numeric_limits<double>::infinity());
  var y2 = std::numeric_limits<double>::infinity();
  var lp2 = stan::math::skew_normal_lcdf(y2, 0.0, 1.0, 0.0);
  EXPECT_EQ(lp2.val(), 0.0);
  lp2.grad();
  EXPECT_EQ(y2.adj(), 0.0);
}
