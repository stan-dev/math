#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <limits>
#include <vector>

class SkewDoubleExponentialTestRig : public VectorRNGTestRig {
 public:
  SkewDoubleExponentialTestRig()
      : VectorRNGTestRig(10000, 10, {-2.5, -1.7, -0.1, 0.1, 2.0},
                         {-3, -2, -1, 0, 2, 6}, {}, {}, {0.1, 1.0, 2.5, 4.0},
                         {1, 2, 3, 4}, {-2.7, -1.5, -0.5, 0.0}, {-3, -2, -1, 0},
                         {0.5}, {0}, {-0.1, 1.1}, {}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& mu, const T2& sigma, const T3& tau,
                        T_rng& rng) const {
    return stan::math::skew_double_exponential_rng(mu, sigma, tau, rng);
  }
};

double icdf(double z, double mu, double sigma, double tau) {
  if (z < tau) {
    return log(z / tau) * sigma / (2.0 * (1.0 - tau)) + mu;
  } else {
    return log((1.0 - z) / (1.0 - tau)) * (-sigma) / (2.0 * tau) + mu;
  }
}

TEST(ProbDistributionsSkewedDoubleExponential, errorCheck) {
  check_dist_throws_all_types(SkewDoubleExponentialTestRig());
}

TEST(ProbDistributionsSkewedDoubleExponential, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::skew_double_exponential_rng(10.0, 2.0, .1, rng));

  EXPECT_THROW(stan::math::skew_double_exponential_rng(10.0, 2.0, -1.0, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::skew_double_exponential_rng(
                   10, 2, stan::math::positive_infinity(), rng),
               std::domain_error);
}

TEST(ProbDistributionsSkewedDoubleExponential, test_sampling_icdf) {
  for (double p : {0.0, 0.1, 0.2, 0.5, 0.7, 0.9, 0.99}) {
    for (double mu : {-1.11, 0.13, 1.2, 4.67}) {
      for (double sigma : {0.11, 1.33}) {
        for (double tau : {0.1, 0.4, 0.77, 0.89}) {
          double x = icdf(p, mu, sigma, tau);
          EXPECT_FLOAT_EQ(
              stan::math::skew_double_exponential_cdf(x, mu, sigma, tau), p);
        }
      }
    }
  }
}

TEST(ProbDistributionsSkewedDoubleExponential, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;

  std::vector<double> samples;
  for (int i = 0; i < N; ++i) {
    samples.push_back(
        stan::math::skew_double_exponential_rng(2.0, 1.0, 0.25, rng));
  }
  std::vector<double> quantiles = {-0.145917216578801,
                                   0.316180903794496,
                                   0.586490975866606,
                                   0.778279024167793,
                                   0.9270413917106,
                                   1.0485890962399,
                                   1.15135621612474,
                                   1.24037714454109,
                                   1.31889916831201,
                                   1.3891395120839,
                                   1.45267963195345,
                                   1.5106872166132,
                                   1.56404902172889,
                                   1.61345433649804,
                                   1.65944958415601,
                                   1.70247526491439,
                                   1.74289167945868,
                                   1.78099728868531,
                                   1.81704210286549,
                                   1.85123763245719,
                                   1.88376440857015,
                                   1.91477775232674,
                                   1.9444122607073,
                                   1.9727853369865,
                                   2,
                                   2.02684604066428,
                                   2.05405734477584,
                                   2.08164398904051,
                                   2.10961647298999,
                                   2.1379857429739,
                                   2.1667632178781,
                                   2.19596081672041,
                                   2.22559098829069,
                                   2.25566674301977,
                                   2.28620168728135,
                                   2.31721006035328,
                                   2.34870677428956,
                                   2.38070745698244,
                                   2.413228498726,
                                   2.44628710262842,
                                   2.47990133926118,
                                   2.51409020597978,
                                   2.54887369140352,
                                   2.58427284560232,
                                   2.62030985660768,
                                   2.65700813394407,
                                   2.69439239996838,
                                   2.73248878990977,
                                   2.77132496162397,
                                   2.81093021621633,
                                   2.85133563085137,
                                   2.89257420525684,
                                   2.9346810236525,
                                   2.97769343409443,
                                   3.02165124753198,
                                   3.0665969592361,
                                   3.1125759956855,
                                   3.15963699050588,
                                   3.20783209366401,
                                   3.25721731884475,
                                   3.30785293481333,
                                   3.35980390761985,
                                   3.41314040178417,
                                   3.4679383501604,
                                   3.52428010409379,
                                   3.5822551778403,
                                   3.64196110413966,
                                   3.70350442147317,
                                   3.76700181810233,
                                   3.83258146374831,
                                   3.90038456709967,
                                   3.97056720672221,
                                   4.04330249506396,
                                   4.11878315102966,
                                   4.19722457733622,
                                   4.27886856637673,
                                   4.36398779521432,
                                   4.45289132035599,
                                   4.54593135162578,
                                   4.64351167996464,
                                   4.74609826873974,
                                   4.85423271128029,
                                   4.96854953896019,
                                   5.08979878259306,
                                   5.2188758248682,
                                   5.3568615678421,
                                   5.50507751214955,
                                   5.66516292749662,
                                   5.83918568147588,
                                   6.02980604108453,
                                   6.24052707240018,
                                   6.47609314371295,
                                   6.743155928962,
                                   7.05145728861651,
                                   7.41610040220442,
                                   7.86238750483284,
                                   8.4377516497364,
                                   9.24868186595273,
                                   10.6349762270726};

  quantiles.push_back(std::numeric_limits<double>::max());

  // Assert that they match
  assert_matches_quantiles(samples, quantiles, 1e-6);
}
