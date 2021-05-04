#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/VectorIntRNGTestRig.hpp>
#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <limits>
#include <vector>

class PoissonLogTestRig : public VectorIntRNGTestRig {
 public:
  PoissonLogTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1, 2, 3, 4, 5, 6},
                            {-0.5, 0.0, 0.1, 1.7}, {-2, 0, 1, 2}, {}, {}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& alpha, const T2&, const T3&,
                        T_rng& rng) const {
    return stan::math::poisson_log_rng(alpha, rng);
  }

  template <typename T1>
  double pmf(int y, T1 alpha, double, double) const {
    return std::exp(stan::math::poisson_log_lpmf(y, alpha));
  }
};

TEST(ProbDistributionsPoissonLog, errorCheck) {
  check_dist_throws_all_types(PoissonLogTestRig());
}

TEST(ProbDistributionsPoissonLog, distributionCheck) {
  check_counts_real(PoissonLogTestRig());
}
