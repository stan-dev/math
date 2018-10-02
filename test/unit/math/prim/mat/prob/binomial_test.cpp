#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/mat/prob/VectorIntRNGTestRig.hpp>
#include <limits>
#include <vector>

class BinomialTestRig : public VectorIntRNGTestRig {
 public:
  BinomialTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1, 2, 3, 4, 5, 6}, {}, {0, 1, 3, 8},
                            {}, {-1, -5, -7}, {0.0, 0.1, 0.7, 0.99}, {0},
                            {-0.1, 1.2}, {-1, 2}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& N, const T2& theta, const T3&,
                        T_rng& rng) const {
    return stan::math::binomial_rng(N, theta, rng);
  }

  template <typename T1>
  double pmf(int y, T1 N, double theta, double) const {
    if (y <= N) {
      return std::exp(stan::math::binomial_lpmf(y, N, theta));
    } else {
      return 0.0;
    }
  }
};

TEST(ProbDistributionsBinomial, errorCheck) {
  check_dist_throws_int_first_argument(BinomialTestRig());
}

TEST(ProbDistributionsBinomial, distributionCheck) {
  check_counts_int_real(BinomialTestRig());
}
