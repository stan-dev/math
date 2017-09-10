#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <limits>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;

struct gumbel_rng_wrapper {
  template<typename T1, typename T2, typename T3, typename T_rng>
  auto operator()(const T1& mu, const T2& sigma, const T3& unused, T_rng& rng) const {
    return stan::math::gumbel_rng(mu, sigma, rng);
  }
};

std::vector<double> build_quantiles(int N, double mu, double sigma,
                                    double unused) {
  std::vector<double> quantiles;
  int K = boost::math::round(2 * std::pow(N, 0.4));
  boost::math::extreme_value_distribution<> dist(mu, sigma);

  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(quantile(dist, frac));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  return quantiles;
}

TEST(ProbDistributionsGumbel, errorCheck) {
  check_dist_throws_all_types(gumbel_rng_wrapper{},
                              {-2.5, -1.7, -0.1, 0.0, 2.0, 5.8}, {},
                              {0.1, 1.0, 2.5, 4.0}, {-2.7, -1.5, -0.5, 0.0},
                              {}, {});
}

TEST(ProbDistributionsGumbel, chiSquareGoodnessFitTest) {
  int N = 10000;
  int M = 10;

  check_quantiles_all_types(N, M, gumbel_rng_wrapper{}, build_quantiles,
                            {-2.5, -1.7, -0.1, 0.0, 2.0, 5.8},
                            {0.1, 1.0, 2.5, 4.0},
                            {});
}
