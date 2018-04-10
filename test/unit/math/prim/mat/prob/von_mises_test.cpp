#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <limits>
#include <vector>

class VonMisesTestRig : public VectorRNGTestRig {
 public:
  VonMisesTestRig()
      : VectorRNGTestRig(
            10000,
            10,
            {-2.5, -1.7, -0.1, 0.0, 2.0, 5.8},
            {-3, -2, -1, 0, 2, 6},
            {}, {}, {0.1, 1.0, 2.5, 4.0}, {1, 2, 3, 4}, {-2.7, -1.5, -0.5, 0.0},
            {-3, -2, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& mu, const T2& kappa, const T3& unused,
                        T_rng& rng) const {
    return stan::math::von_mises_rng(mu, kappa, rng);
  }
};

TEST(ProbDistributionsVonMises, errorCheck) {
  check_dist_throws_all_types(VonMisesTestRig());
}
