#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <limits>
#include <vector>

class UniformTestRig : public VectorRealRNGTestRig {
 public:
  UniformTestRig(std::vector<double> good_p1, std::vector<int> good_p1_int,
                 std::vector<double> bad_p1, std::vector<int> bad_p1_int,
                 std::vector<double> good_p2, std::vector<int> good_p2_int,
                 std::vector<double> bad_p2, std::vector<int> bad_p2_int)
      : VectorRealRNGTestRig(10000, 10, good_p1, good_p1_int, bad_p1,
                             bad_p1_int, good_p2, good_p2_int, bad_p2,
                             bad_p2_int) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& alpha, const T2& beta, const T3& unused,
                        T_rng& rng) const {
    return stan::math::uniform_rng(alpha, beta, rng);
  }

  std::vector<double> generate_quantiles(double alpha, double beta,
                                         double unused) const {
    std::vector<double> quantiles;
    double K = stan::math::round(2 * std::pow(N_, 0.4));
    boost::math::uniform_distribution<> dist(alpha, beta);

    for (int i = 1; i < K; ++i) {
      double frac = i / K;
      quantiles.push_back(quantile(dist, frac));
    }
    quantiles.push_back(std::numeric_limits<double>::max());

    return quantiles;
  }
};

TEST(ProbDistributionsUniform, errorCheck) {
  // beta must be greater than alpha, so we have to be careful about how the
  // tests are initialized here (the VectorRNGTestRig class assumes no order)

  check_dist_throws_all_types(
      UniformTestRig({-2.5, -1.7, -0.1, 0.0, 0.5}, {-3, -2, -1, 0, 1}, {}, {},
                     {1.1, 2.2, 3.8}, {2, 3, 6}, {}, {}));

  check_dist_throws_all_types(
      UniformTestRig({-7.5, -6.7, -5.1, -4.0}, {-7, -6, -5, -4}, {}, {},
                     {-3.0, -2.2, 0.0, 1.0}, {-3, -2, 0, 1, 2}, {}, {}));
}

TEST(ProbDistributionsUniform, distributionTest) {
  check_quantiles_real_real(UniformTestRig({-1.7, -0.1, 0.0, 0.5},
                                           {-2, -1, 0, 1}, {}, {}, {1.1, 3.8},
                                           {2, 6}, {}, {}));

  check_quantiles_real_real(UniformTestRig({-7.5, -5.1, -4.0}, {-7, -6, -4}, {},
                                           {}, {-2.2, 0.0, 1.0}, {-3, 0, 1, 2},
                                           {}, {}));
}
