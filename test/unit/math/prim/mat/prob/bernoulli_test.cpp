#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/mat/prob/VectorIntRNGTestRig.hpp>
#include <limits>
#include <vector>

class BernoulliTestRig : public VectorIntRNGTestRig {
 public:
  BernoulliTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1}, {0.0, 0.1, 0.2, 0.7, 1.0},
                            {0, 1}, {-2.0, -0.5, 1.1, 2.0}, {-2, -1, 2}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& theta, const T2&, const T3&,
                        T_rng& rng) const {
    return stan::math::bernoulli_rng(theta, rng);
  }

  template <typename T1>
  double pmf(int y, T1 theta, double, double) const {
    return std::exp(stan::math::bernoulli_lpmf(y, theta));
  }
};

TEST(ProbDistributionsBernoulli, errorCheck) {
  check_dist_throws_all_types(BernoulliTestRig());
}

TEST(ProbDistributionsBernoulli, distributionCheck) {
  check_counts_real(BernoulliTestRig());
}
