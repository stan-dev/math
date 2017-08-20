#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>

using Eigen::Dynamic;
using Eigen::Matrix;

struct dist_wrapper {
  template<typename T1, typename T2, typename T3>
  auto operator()(const T1& a1, const T2& a2, T3& a3) const {
    return stan::math::normal_rng(a1, a2, a3);
  }
};

std::vector<double> build_quantiles(int N, double mu, double sigma) {
  std::vector<double> quantiles;
  int K = boost::math::round(2 * std::pow(N, 0.4));
  boost::math::normal_distribution<> dist(mu, sigma);

  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(quantile(dist, frac));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  return quantiles;
};

TEST(ProbDistributionsNormal, error_check) {
  twoParamDistCheckThrowsAllTypes(dist_wrapper{},
                                  Constraint::None, Constraint::Positive);
}

TEST(ProbDistributionsNormal, chiSquareGoodnessFitTest) {
  int N = 10000;
  int M = 10;

  chiSquareTestAllTypes(N, M, dist_wrapper{}, build_quantiles,
                        Constraint::None, Constraint::Positive);
}
