#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/mat/prob/VectorIntRNGTestRig.hpp>
#include <limits>
#include <vector>

class BernoulliLogitTestRig : public VectorIntRNGTestRig {
 public:
  BernoulliLogitTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1},
                            {-5.7, -1.0, 0.0, 0.2, 1.0, 10.0}, {-3, -2, 0, 1},
                            {}, {}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& t, const T2&, const T3&, T_rng& rng) const {
    return stan::math::bernoulli_logit_rng(t, rng);
  }

  template <typename T1>
  double pmf(int y, T1 t, double, double) const {
    return std::exp(stan::math::bernoulli_logit_lpmf(y, t));
  }
};

TEST(ProbDistributionsBernoulliLogit, errorCheck) {
  check_dist_throws_all_types(BernoulliLogitTestRig());
}

TEST(ProbDistributionsBernoulliLogit, distributionCheck) {
  check_counts_real(BernoulliLogitTestRig());
}
