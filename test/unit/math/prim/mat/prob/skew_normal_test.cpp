#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <limits>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;

struct skew_normal_rng_wrapper {
  template<typename T1, typename T2, typename T3, typename T_rng>
  auto operator()(const T1& mu, const T2& sigma, const T3& alpha,
                  T_rng& rng) const {
    return stan::math::skew_normal_rng(mu, sigma, alpha, rng);
  }
};

std::vector<double> build_quantiles(int N, double mu, double sigma,
                                    double alpha) {
  std::vector<double> quantiles;
  int K = boost::math::round(2 * std::pow(N, 0.4));
  boost::math::skew_normal_distribution<> dist(mu, sigma, alpha);

  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(quantile(dist, frac));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  return quantiles;
}

TEST(ProbDistributionsSkewNormal, errorCheck) {
  check_dist_throws_all_types(skew_normal_rng_wrapper{},
                              {-2.5, -1.7, -0.1, 0.0, 2.0, 5.8}, {},
                              {0.1, 1.0, 2.5, 4.0}, {-2.7, -1.5, -0.5, 0.0},
                              {-2.5, -1.7, -1.3, 0.0, 1.0, 4.8}, {});
}

TEST(ProbDistributionsSkewNormal, chiSquareGoodnessFitTest) {
  int N = 10000;
  int M = 10;

  check_quantiles_all_types(N, M, skew_normal_rng_wrapper{}, build_quantiles,
                              {-2.5, -1.7, -0.1, 0.0, 2.0, 5.8},
                              {0.1, 1.0, 2.5, 4.0},
                              {-2.5, -1.7, -1.3, 0.0, 1.0, 4.8});
}

