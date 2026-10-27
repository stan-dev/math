#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Values and gradients of beta_binomial_{cdf,lcdf,lccdf} at points that
// cover the routes of internal::beta_binomial_log_cdfs(): the lower or the
// upper tail summed first, the other one as its complement or summed too,
// then the first one as the complement of a small second one, the unimodal
// and the U-shaped pmf, an early stop, tails near 1e-276, and 200 000
// trials.
//
// References: mpmath at 50 digits; both tails are direct sums of the pmf,
// with the derivatives as sums of pmf * d log pmf; the larger tail is then
// taken as 1 minus the smaller one, so that its log keeps its digits. The
// shapes are written in hex so that they are exact. The partials are
// compared as x * d/dx, the gradient in log(x) that a sampler sees.

namespace beta_binomial_log_cdfs_test_internal {

struct TestCase {
  int n;
  int N;
  double alpha;
  double beta;
  double lcdf;
  double lccdf;
  std::vector<double> dlcdf;   // d/d(alpha, beta)
  std::vector<double> dlccdf;  // d/d(alpha, beta)
};

// clang-format off
// NOLINTBEGIN(whitespace/line_length)
const std::vector<TestCase> test_cases = {
    {2, 10, 0x1.8000000000000p+1, 0x1.0000000000000p+2, -1.4198170531585341, -0.27675300191959068, {-0.51879717504717504, 0.30789835164835166}, {0.16541359204402684, -0.098170488931358491}},  // lower tail first, upper as its complement
    {8, 10, 0x1.8000000000000p+1, 0x1.0000000000000p+2, -0.036367644170874847, -3.3322045101752038, {-0.022532560032560031, 0.027348937071159294}, {0.60837912087912083, -0.73842130092130087}},  // upper tail first, lower as its complement
    {8, 10, 0x1.0000000000000p-1, 0x1.e000000000000p+4, -7.8114337516020798e-09, -18.667677314898757, {-3.0283500032506201e-08, 1.994662475956812e-09}, {3.8768170962233981, -0.25535164626559514}},  // small upper tail
    {1, 10, 0x1.e000000000000p+4, 0x1.0000000000000p-1, -18.667677314898757, -7.8114337516020798e-09, {-0.25535164626559514, 3.8768170962233981}, {1.994662475956812e-09, -3.0283500032506201e-08}},  // small lower tail
    {58, 117, 0x1.0000000000000p+1, 0x1.0000000000000p+1, -0.69314718055994529, -0.69314718055994529, {-0.43689703406552266, 0.43689703406552266}, {0.43689703406552266, -0.43689703406552266}},  // lower tail above 1/2: both tails summed
    {58, 117, 0x1.4000000000000p+3, 0x1.3880000000000p+13, -4.4746105258856488e-112, -256.39111110071553, {-8.7861915440757089e-112, 2.6123844465240352e-114}, {1.9635656540938118, -0.0058382387280666677}},  // both tails summed, lcdf as the complement of the upper tail
    {59, 117, 0x1.3880000000000p+13, 0x1.4000000000000p+3, -251.39436575466911, -6.6193321069869591e-110, {-0.005738828270701637, 1.9488631096182969}, {3.7987210228739794e-112, -1.2900172153618836e-109}},  // both tails summed, lccdf as the complement of the lower tail
    {60, 117, 0x1.3333333333333p-2, 0x1.999999999999ap-3, -0.91507251589898242, -0.51163859300372261, {-2.1340239197768334, 3.1499047202585775}, {1.4255752818494947, -2.1042061749013872}},  // U-shaped pmf: no early stop
    {5, 117, 0x1.f400000000000p+9, 0x1.7700000000000p+11, -19.695873917018684, -2.7937699528797345e-09, {-0.023991525922721486, 0.0078730012808713602}, {6.7026804340263269e-11, -2.1995354448207051e-11}},  // large shapes: early stop of the lower sum
    {3, 10, 0x1.2a4da7abab141p+8, 0x1.258191c361692p+7, -3.933350026668649, -0.019772163107297, {-0.012878429074005715, 0.025788239904711231}, {0.00025716841006185448, -0.00051496348013239945}},  // the 3F2 series once gave NaN here
    {0, 117, 0x1.3880000000000p+13, 0x1.0000000000000p+0, -634.8223474979892, -1.9959860090434558e-276, {-0.011631506014216147, 5.3320255103901211}, {2.3216323268480239e-278, -1.0642648318601472e-275}},  // cdf about 2e-276
    {116, 117, 0x1.0000000000000p+0, 0x1.3880000000000p+13, -1.9959860090434558e-276, -634.8223474979892, {-1.0642648318601472e-275, 2.3216323268480239e-278}, {5.3320255103901211, -0.011631506014216147}},  // ccdf about 2e-276
    {0, 1, 0x1.0000000000000p+1, 0x1.8000000000000p+1, -0.51082562376599072, -0.91629073187415511, {-0.20000000000000001, 0.13333333333333333}, {0.29999999999999999, -0.20000000000000001}},  // N = 1
    {1000, 10000, 0x1.9000000000000p+6, 0x1.e000000000000p+4, -162.91049486808714, -1.773662707038273e-71, {-1.9664909412961633, 1.3654243090365596}, {3.4878916463055947e-71, -2.421802176221648e-71}},  // lccdf gradient about 1e-71
    {50000, 200000, 0x1.0000000000000p-1, 0x1.0000000000000p-1, -1.0986053971070815, -0.40546855390648867, {-1.9383781898621288, 1.2922610759880457}, {0.96919911383880364, -0.64613721731210594}},  // 200 000 trials, U-shaped
};
// NOLINTEND
// clang-format on

}  // namespace beta_binomial_log_cdfs_test_internal

TEST_F(AgradRev, ProbDistributionsBetaBinomialCdfs_values_and_gradients) {
  using beta_binomial_log_cdfs_test_internal::test_cases;
  using stan::math::var;
  // the values relative to the reference, so that an lcdf or lccdf near 0
  // keeps its digits; the gradients relative to max(1, |reference|)
  constexpr double tol_value = 1e-11;
  constexpr double tol_grad = 1e-10;
  for (const auto& c : test_cases) {
    const std::vector<double> x{c.alpha, c.beta};
    for (int which = 0; which < 3; ++which) {
      var alpha = c.alpha;
      var beta = c.beta;
      var lp;
      if (which == 0) {
        lp = stan::math::beta_binomial_cdf(c.n, c.N, alpha, beta);
      } else if (which == 1) {
        lp = stan::math::beta_binomial_lcdf(c.n, c.N, alpha, beta);
      } else {
        lp = stan::math::beta_binomial_lccdf(c.n, c.N, alpha, beta);
      }
      lp.grad();
      const std::vector<double> grads{alpha.adj(), beta.adj()};
      // the partials of the cdf are cdf * d lcdf
      const double cdf = std::exp(c.lcdf);
      const double ref = (which == 0) ? cdf : (which == 1) ? c.lcdf : c.lccdf;
      const std::vector<double>& dref = (which == 2) ? c.dlccdf : c.dlcdf;
      const double scale = (which == 0) ? cdf : 1.0;
      EXPECT_NEAR(lp.val(), ref, tol_value * std::fabs(ref))
          << "function " << which << ", n = " << c.n << ", N = " << c.N
          << ", alpha = " << c.alpha << ", beta = " << c.beta;
      for (int i = 0; i < 2; ++i) {
        const double log_scale_ref = x[i] * dref[i];
        EXPECT_NEAR(x[i] * grads[i], scale * log_scale_ref,
                    tol_grad * scale * std::max(1.0, std::fabs(log_scale_ref)))
            << "function " << which << ", gradient " << i << ", n = " << c.n
            << ", N = " << c.N << ", alpha = " << c.alpha
            << ", beta = " << c.beta;
      }
      stan::math::recover_memory();
    }
  }
}

TEST_F(AgradRev, ProbDistributionsBetaBinomialCdfs_vectors_and_edges) {
  using stan::math::beta_binomial_cdf;
  using stan::math::beta_binomial_lccdf;
  using stan::math::beta_binomial_lcdf;
  using stan::math::var;
  const std::vector<int> ns{3, 40, 0};
  const std::vector<int> Ns{10, 117, 5};
  const std::vector<double> alphas{2.5, 30.0, 0.5};
  const double beta = 7.0;

  // the lcdf and lccdf of vectors are sums, the cdf is a product
  double lcdf_sum = 0;
  double lccdf_sum = 0;
  for (size_t i = 0; i < ns.size(); ++i) {
    lcdf_sum += beta_binomial_lcdf(ns[i], Ns[i], alphas[i], beta);
    lccdf_sum += beta_binomial_lccdf(ns[i], Ns[i], alphas[i], beta);
  }
  EXPECT_NEAR(beta_binomial_lcdf(ns, Ns, alphas, beta), lcdf_sum, 1e-14);
  EXPECT_NEAR(beta_binomial_lccdf(ns, Ns, alphas, beta), lccdf_sum, 1e-14);
  EXPECT_NEAR(beta_binomial_cdf(ns, Ns, alphas, beta), std::exp(lcdf_sum),
              1e-14 * std::exp(lcdf_sum));

  // the gradient of a vector argument holds the partials of its elements
  std::vector<var> av(alphas.begin(), alphas.end());
  var lp = beta_binomial_lccdf(ns, Ns, av, beta);
  lp.grad();
  std::vector<double> adj_vector;
  for (const var& v : av) {
    adj_vector.push_back(v.adj());
  }
  for (size_t i = 0; i < ns.size(); ++i) {
    var ai = alphas[i];
    var lpi = beta_binomial_lccdf(ns[i], Ns[i], ai, beta);
    stan::math::set_zero_all_adjoints();
    lpi.grad();
    EXPECT_NEAR(adj_vector[i], ai.adj(),
                1e-14 * std::max(1.0, std::fabs(ai.adj())));
  }
  stan::math::recover_memory();

  // n < 0 is a term of log(1) = 0 in the lccdf and makes the lcdf -inf;
  // n >= N is a term of 0 in the lcdf and makes the lccdf -inf
  const double inf = std::numeric_limits<double>::infinity();
  EXPECT_DOUBLE_EQ(beta_binomial_lccdf(std::vector<int>{-1, 7},
                                       std::vector<int>{20, 20}, 2.5, 4.0),
                   beta_binomial_lccdf(7, 20, 2.5, 4.0));
  EXPECT_EQ(beta_binomial_lcdf(std::vector<int>{-1, 7},
                               std::vector<int>{20, 20}, 2.5, 4.0),
            -inf);
  EXPECT_DOUBLE_EQ(beta_binomial_lcdf(std::vector<int>{20, 7},
                                      std::vector<int>{20, 20}, 2.5, 4.0),
                   beta_binomial_lcdf(7, 20, 2.5, 4.0));
  EXPECT_EQ(beta_binomial_lccdf(std::vector<int>{20, 7},
                                std::vector<int>{20, 20}, 2.5, 4.0),
            -inf);
}
