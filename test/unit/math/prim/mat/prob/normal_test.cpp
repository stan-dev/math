#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <limits>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;

/*
 * This needs to be a functor instead of a function because we need to pass it
 * as an argument to another function (and I'm not sure it's possible to easily
 * pass a templated function as an argument)
 */
struct normal_rng_wrapper {
  template <typename T1, typename T2, typename T3>
  auto operator()(const T1& mean, const T2& sd, T3& rng) const {
    return stan::math::normal_rng(mean, sd, rng);
  }
};

/*
 * This function builds the quantiles that we will supply to
 * assert_matches_quantiles to test the normal_rng
 */
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
}

TEST(ProbDistributionsNormal, errorCheck) {
  /*
   * This test verifies that normal_rng throws errors in the right places. Test
   * inputs are picked from the four input initializer lists which are (in
   * order): valid values for p1, invalid values for p1,
   *         valid values for p2, and invalid values for p2.
   *
   * It does so for all possible combinations of calling arguments.
   */
  check_dist_throws_all_types(normal_rng_wrapper{},
                              {-2.5, -1.7, -0.1, 0.0, 2.0, 5.8}, {},
                              {0.1, 1.0, 2.5, 4.0}, {-2.7, -1.5, -0.5, 0.0});
}

TEST(ProbDistributionsNormal, chiSquareGoodnessFitTest) {
  int N = 10000;
  int M = 10;

  /*
   * This test checks that the normal_rng is actually generating numbers from
   * the correct distributions. Test inputs are picked from the two input
   * initializer lists which are (in order): valid values for p1, and valid
   * values for p2.
   *
   * It does so for all possible combinations of calling arguments.
   */
  check_quantiles_all_types(N, M, normal_rng_wrapper{}, build_quantiles,
                            {-2.5, -1.7, -0.1, 0.0, 2.0, 5.8},
                            {0.1, 1.0, 2.5, 4.0});
}
