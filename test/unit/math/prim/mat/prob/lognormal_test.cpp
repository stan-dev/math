#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <limits>
#include <vector>

class LogNormalTestRig : public VectorRealRNGTestRig {
 public:
  LogNormalTestRig()
      : VectorRealRNGTestRig(10000, 10, {-2.5, -1.7, -0.1, 0.0, 2.0, 5.8},
                             {-3, -2, -1, 0, 2, 6}, {}, {},
                             {0.1, 1.0, 2.5, 4.0}, {1, 2, 3, 4},
                             {-2.7, -1.5, -0.5, 0.0}, {-3, -2, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& mean, const T2& sd, const T3& unused,
                        T_rng& rng) const {
    return stan::math::lognormal_rng(mean, sd, rng);
  }

  std::vector<double> generate_quantiles(double mu, double sigma,
                                         double unused) const {
    std::vector<double> quantiles;
    double K = stan::math::round(2 * std::pow(N_, 0.4));
    boost::math::lognormal_distribution<> dist(mu, sigma);

    for (int i = 1; i < K; ++i) {
      double frac = i / K;
      quantiles.push_back(quantile(dist, frac));
    }
    quantiles.push_back(std::numeric_limits<double>::max());

    return quantiles;
  }
};

TEST(ProbDistributionsLogNormal, errorCheck) {
  check_dist_throws_all_types(LogNormalTestRig());
}

TEST(ProbDistributionsLogNormal, chiSquareGoodnessFitTest) {
  check_quantiles_real_real(LogNormalTestRig());
}
