#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Values and gradients of beta_neg_binomial_{cdf,lcdf,lccdf} at points that
// cover the routes of internal::beta_neg_binomial_log_cdfs(): the lower sum
// from n down, the complement of either tail, the short and the long tail
// sum, and the Thomae series with a non-integer and with an integer x. The
// points include small alpha, where the terms of the upper tail fall only
// like k^(-1 - alpha); shapes of 1e6; a cdf of 2.6e-35; and n = 20000.
//
// References: mpmath at 50 digits or more. The cdf is the finite sum of the
// pmf; the ccdf is 1 - cdf computed with extra digits, or the sum of the pmf
// over k > n; the gradients are the sums of pmf * d log pmf, checked against
// mp.diff. The parameters are written in hex so that they are exact. The
// partials are compared as x * d/dx, the gradient in log(x) that a sampler
// sees.

namespace beta_neg_binomial_cdfs_test_internal {

struct TestCase {
  int n;
  double r;
  double alpha;
  double beta;
  double lcdf;
  double lccdf;
  std::vector<double> dlcdf;   // d/d(r, alpha, beta)
  std::vector<double> dlccdf;  // d/d(r, alpha, beta)
};

// clang-format off
// NOLINTBEGIN(whitespace/line_length)
const std::vector<TestCase> test_cases = {
    {2, 0x1.8000000000000p+1, 0x1.0000000000000p+2, 0x1.8000000000000p+1, -0.49721022174799306, -0.93706785981591312, {-0.19621186445029862, 0.15355088762562072, -0.19621186445029862}, {0.30461620945046358, -0.23838563217016256, 0.30461620945046358}},  // ordinary point
    {20, 0x1.8000000000000p+1, 0x1.0000000000000p-1, 0x1.8000000000000p+1, -0.87140397613822818, -0.54191005068750497, {-0.1891559327062384, 1.6158672896523629, -0.1891559327062384}, {0.13605752717370073, -1.1622733927801345, 0.13605752717370073}},  // small alpha: the 3F2 value was wrong by 7e-3
    {20, 0x1.8000000000000p+1, 0x1.0000000000000p-2, 0x1.8000000000000p+1, -1.4760580874852567, -0.25946629095556611, {-0.23559870423820939, 3.6843570031945654, -0.23559870423820939}, {0.069793309685380625, -1.0914468742386048, 0.069793309685380625}},  // upper tail first, its complement for the cdf
    {20000, 0x1.8000000000000p+1, 0x1.0000000000000p-2, 0x1.8000000000000p+1, -0.16335403587918723, -1.8924008435000104, {-0.016711519065541223, 1.3561824508346239, -0.016711519065541223}, {0.094174086515496705, -7.6424676269584797, 0.094174086515496705}},  // heavy tail, large n
    {20000, 0x1.0000000000000p-1, 0x1.0000000000000p+0, 0x1.0000000000000p-1, -1.2499375029946549e-05, -11.289838162191222, {-2.499843759374463e-05, 0.00012816151869144557, -2.499843759374463e-05}, {1.999962502265497, -10.25337006183252, 1.999962502265497}},  // Thomae series, x = 0.5
    {0, 0x1.0000000000000p-1, 0x1.0000000000000p-2, 0x1.0000000000000p-1, -0.78318878541367354, -0.61054758643660245, {-0.85840734641020677, 2.2831853071795867, -0.85840734641020677}, {0.72229782263411269, -1.921162234855734, 0.72229782263411269}},  // n = 0, small shapes
    {0, 0x1.e000000000000p+4, 0x1.4000000000000p+3, 0x1.2c00000000000p+8, -79.627954609150891, -2.6182957935133381e-35, {-2.1511469344641077, 1.3320590035711872, -0.09251578139693474}, {5.6323389697364848e-35, -3.4877244857620081e-35, 2.4223368126519377e-36}},  // cdf 2.6e-35: 1 - ccdf gave a NaN gradient
    {0, 0x1.e848000000000p+19, 0x1.0000000000000p-2, 0x1.e848000000000p+19, -1386298.010794742, 0, {-0.69314730555993753, 17.349816535780612, -0.69314730555993753}, {0, 0, 0}},  // shapes of 1e6: the 3F2 value overflowed
    {2, 0x1.e848000000000p+19, 0x1.0000000000000p+2, 0x1.e848000000000p+19, -1386223.7540788064, 0, {-0.69314443056525776, 11.866249958966133, -0.69314443056525776}, {0, 0, 0}},  // shapes of 1e6: lcdf -1.4e6, lower sum from n down
    {200, 0x1.8000000000000p+1, 0x1.3880000000000p+13, 0x1.e000000000000p+4, 0, -890.84885316612292, {0, 0, 0}, {4.370163163488674, -0.019838211951244258, 2.0357763248986651}},  // tail sum: lccdf -891
    {20000, 0x1.8000000000000p+1, 0x1.9000000000000p+6, 0x1.d4c0000000000p+14, -5.8901075871069618e-20, -44.27842759632194, {-1.6515719040741391e-19, 2.8870435083041056e-20, -7.6037348214479634e-23}, {2.8039757842272959, -0.49015123503408398, 0.0012909330957030416}},  // Thomae series, integer x = 3
    {20000, 0x1.e000000000000p+4, 0x1.9000000000000p+6, 0x1.d4c0000000000p+14, -1.9902072448370877e-05, -10.824696639541852, {-1.1663144244205176e-05, 5.2405299430920191e-06, -1.5404356555183488e-08}, {0.58602078826087212, -0.26331317043269375, 0.00077399995936819138}},  // long tail sum (5.7e4 terms)
    {2000, 0x1.e000000000000p+4, 0x1.0000000000000p+1, 0x1.2c00000000000p+8, -2.5815723689661887, -0.078669849505074954, {-0.11166547403255754, 1.1893396752518539, -0.010649489062478583}, {0.0091394941869035913, -0.09734399233418238, 0.00087162969774917404}},  // the 3F2 value was wrong by 5e-4
    {200, 0x1.2c00000000000p+8, 0x1.4000000000000p+3, 0x1.7700000000000p+11, -649.56328188146324, 0, {-1.9181174130266818, 3.3297281096790732, -0.089172183391557094}, {1.5174913349159807e-282, -2.6342670786201108e-282, 7.0547305755752653e-284}},  // mid shapes, alpha 10
    {20, 0x1.d4c0000000000p+14, 0x1.3880000000000p+13, 0x1.0000000000000p-1, -0.00054974363964707012, -7.5063333573606696, {-1.0207926638369315e-07, 3.0634751461336257e-07, -0.0020553852035363543}, {0.00018563418088173196, -0.0005571020634753558, 3.7377791021806801}},  // ccdf 5.5e-4: 1 - cdf would lose 3 digits
};
// NOLINTEND
// clang-format on

}  // namespace beta_neg_binomial_cdfs_test_internal

TEST_F(AgradRev, ProbDistributionsBetaNegBinomialCdfs_values_and_gradients) {
  using beta_neg_binomial_cdfs_test_internal::test_cases;
  using stan::math::var;
  // tolerances times max(1, |reference|); relative for the cdf
  constexpr double tol_value = 1e-11;
  constexpr double tol_grad = 1e-10;
  for (const auto& c : test_cases) {
    const std::vector<double> x{c.r, c.alpha, c.beta};
    for (int which = 0; which < 3; ++which) {
      var r = c.r;
      var alpha = c.alpha;
      var beta = c.beta;
      var lp;
      if (which == 0) {
        lp = stan::math::beta_neg_binomial_cdf(c.n, r, alpha, beta);
      } else if (which == 1) {
        lp = stan::math::beta_neg_binomial_lcdf(c.n, r, alpha, beta);
      } else {
        lp = stan::math::beta_neg_binomial_lccdf(c.n, r, alpha, beta);
      }
      lp.grad();
      const std::vector<double> grads{r.adj(), alpha.adj(), beta.adj()};
      // the partials of the cdf are cdf * d lcdf
      const double cdf = std::exp(c.lcdf);
      const double ref = (which == 0) ? cdf : (which == 1) ? c.lcdf : c.lccdf;
      const std::vector<double>& dref = (which == 2) ? c.dlccdf : c.dlcdf;
      const double scale = (which == 0) ? cdf : 1.0;
      const double value_tol = (which == 0)
                                   ? tol_value * cdf
                                   : tol_value * std::max(1.0, std::fabs(ref));
      EXPECT_NEAR(lp.val(), ref, value_tol)
          << "function " << which << ", n = " << c.n << ", r = " << c.r
          << ", alpha = " << c.alpha << ", beta = " << c.beta;
      for (int i = 0; i < 3; ++i) {
        const double log_scale_ref = x[i] * dref[i];
        EXPECT_NEAR(x[i] * grads[i], scale * log_scale_ref,
                    tol_grad * scale * std::max(1.0, std::fabs(log_scale_ref)))
            << "function " << which << ", gradient " << i << ", n = " << c.n
            << ", r = " << c.r << ", alpha = " << c.alpha
            << ", beta = " << c.beta;
      }
      stan::math::recover_memory();
    }
  }
}

TEST_F(AgradRev, ProbDistributionsBetaNegBinomialCdfs_vectors_and_largest_int) {
  using stan::math::beta_neg_binomial_cdf;
  using stan::math::beta_neg_binomial_lccdf;
  using stan::math::beta_neg_binomial_lcdf;
  using stan::math::var;
  const std::vector<int> ns{3, 40, 0};
  const std::vector<double> rs{2.5, 30.0, 0.5};
  const double alpha = 3.0;
  const double beta = 7.0;

  // the lcdf and lccdf of vectors are sums, the cdf is a product
  double lcdf_sum = 0;
  double lccdf_sum = 0;
  for (size_t i = 0; i < ns.size(); ++i) {
    lcdf_sum += beta_neg_binomial_lcdf(ns[i], rs[i], alpha, beta);
    lccdf_sum += beta_neg_binomial_lccdf(ns[i], rs[i], alpha, beta);
  }
  EXPECT_NEAR(beta_neg_binomial_lcdf(ns, rs, alpha, beta), lcdf_sum, 1e-14);
  EXPECT_NEAR(beta_neg_binomial_lccdf(ns, rs, alpha, beta), lccdf_sum, 1e-14);
  EXPECT_NEAR(beta_neg_binomial_cdf(ns, rs, alpha, beta), std::exp(lcdf_sum),
              1e-14 * std::exp(lcdf_sum));

  // the gradient of a vector argument holds the partials of its elements
  std::vector<var> rv(rs.begin(), rs.end());
  var lp = beta_neg_binomial_lccdf(ns, rv, alpha, beta);
  lp.grad();
  std::vector<double> adj_vector;
  for (const var& v : rv) {
    adj_vector.push_back(v.adj());
  }
  for (size_t i = 0; i < ns.size(); ++i) {
    var ri = rs[i];
    var lpi = beta_neg_binomial_lccdf(ns[i], ri, alpha, beta);
    stan::math::set_zero_all_adjoints();
    lpi.grad();
    EXPECT_NEAR(adj_vector[i], ri.adj(),
                1e-14 * std::max(1.0, std::fabs(ri.adj())));
  }
  stan::math::recover_memory();

  // the largest int stands for infinity: a cdf factor of 1, an lcdf term of
  // 0, and an lccdf of -inf
  const std::vector<int> n_inf{std::numeric_limits<int>::max(), 5};
  EXPECT_DOUBLE_EQ(beta_neg_binomial_lcdf(n_inf, 3.0, 2.0, 4.0),
                   beta_neg_binomial_lcdf(5, 3.0, 2.0, 4.0));
  EXPECT_DOUBLE_EQ(beta_neg_binomial_cdf(n_inf, 3.0, 2.0, 4.0),
                   beta_neg_binomial_cdf(5, 3.0, 2.0, 4.0));
  EXPECT_EQ(beta_neg_binomial_lccdf(n_inf, 3.0, 2.0, 4.0),
            -std::numeric_limits<double>::infinity());
}
