#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/mat/prob/VectorIntRNGTestRig.hpp>
#include <limits>
#include <vector>

class NegativeBinomialTestRig : public VectorIntRNGTestRig {
 public:
  NegativeBinomialTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1, 2, 3, 4, 5, 6}, {0.1, 1.7, 3.99},
                            {1, 2, 3}, {-2.1, -0.5, 0.0}, {-3, -1, 0},
                            {0.1, 1.1, 4.99}, {1, 2, 3}, {-3.0, -2.0, 0.0},
                            {-3, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& alpha, const T2& beta, const T3&,
                        T_rng& rng) const {
    return stan::math::neg_binomial_rng(alpha, beta, rng);
  }

  template <typename T1>
  double pmf(int y, T1 alpha, double beta, double) const {
    return std::exp(stan::math::neg_binomial_lpmf(y, alpha, beta));
  }
};

TEST(ProbDistributionsNegativeBinomial, errorCheck) {
  check_dist_throws_all_types(NegativeBinomialTestRig());
}

TEST(ProbDistributionsNegativeBinomial, distributionCheck) {
  check_counts_real_real(NegativeBinomialTestRig());
}
