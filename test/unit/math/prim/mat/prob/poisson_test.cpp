#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/mat/prob/VectorIntRNGTestRig.hpp>
#include <limits>
#include <vector>

class PoissonTestRig : public VectorIntRNGTestRig {
 public:
  PoissonTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1, 2, 3, 4, 5, 6}, {0.1, 1.1, 4.99},
                            {1, 2, 3}, {-3.0, -2.0, 0.0}, {-3, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& lambda, const T2&, const T3&,
                        T_rng& rng) const {
    return stan::math::poisson_rng(lambda, rng);
  }

  template <typename T1>
  double pmf(int y, T1 lambda, double, double) const {
    return std::exp(stan::math::poisson_lpmf(y, lambda));
  }
};

TEST(ProbDistributionsPoisson, errorCheck) {
  check_dist_throws_all_types(PoissonTestRig());
}

TEST(ProbDistributionsPoisson, distributionCheck) {
  check_counts_real(PoissonTestRig());
}
