#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/mat/prob/VectorIntRNGTestRig.hpp>
#include <limits>
#include <vector>

class NegativeBinomial2TestRig : public VectorIntRNGTestRig {
 public:
  NegativeBinomial2TestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1, 2, 3, 4, 5, 6}, {0.1, 1.7, 3.99},
                            {1, 2, 3}, {-2.1, -0.5, 0.0}, {-3, -1, 0},
                            {0.1, 1.1, 4.99}, {1, 2, 3}, {-3.0, -2.0, 0.0},
                            {-3, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& mu, const T2& phi, const T3&,
                        T_rng& rng) const {
    return stan::math::neg_binomial_2_rng(mu, phi, rng);
  }

  template <typename T1>
  double pmf(int y, T1 mu, double phi, double) const {
    return std::exp(stan::math::neg_binomial_2_lpmf(y, mu, phi));
  }
};

TEST(ProbDistributionsNegativeBinomial2, errorCheck) {
  check_dist_throws_all_types(NegativeBinomial2TestRig());
}

TEST(ProbDistributionsNegativeBinomial2, distributionCheck) {
  check_counts_real_real(NegativeBinomial2TestRig());
}
