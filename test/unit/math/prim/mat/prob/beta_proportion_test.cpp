#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <limits>
#include <vector>

class BetaProportionTestRig : public VectorRealRNGTestRig {
 public:
  BetaProportionTestRig()
      : VectorRealRNGTestRig(10000, 10, {0.3, 0.4, 0.5, 0.6, 0.7}, {1, 2, 3},
                             {-2.5, -1.7, -0.1, 0.0}, {-3, -2, -1, 0},
                             {0.35, 0.5, 0.9, 1.7, 2.1, 4.1}, {1, 2, 3, 4},
                             {-2.7, -1.5, -0.5, 0.0}, {-3, -2, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& mu, const T2& kappa, const T3&,
                        T_rng& rng) const {
    return stan::math::beta_proportion_rng(mu, kappa, rng);
  }

  std::vector<double> generate_quantiles(double mu, double kappa,
                                         double) const {
    // transform from location and precision parameterization
    // into shape1 (alpha) and shape2 (beta) parameterization
    double alpha = mu * kappa;
    double beta = kappa - alpha;
    std::vector<double> quantiles;
    double K = stan::math::round(2 * std::pow(N_, 0.4));
    boost::math::beta_distribution<> dist(alpha, beta);

    for (int i = 1; i < K; ++i) {
      double frac = i / K;
      quantiles.push_back(quantile(dist, frac));
    }
    quantiles.push_back(std::numeric_limits<double>::max());

    return quantiles;
  }
};

TEST(ProbDistributionsBetaProportion, errorCheck) {
  check_dist_throws_real_first_argument(BetaProportionTestRig());
}

TEST(ProbDistributionsBetaProportion, distributionTest) {
  check_quantiles_real_first_argument(BetaProportionTestRig());
}
