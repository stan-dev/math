#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <vector>

class ParetoType2TestRig : public VectorRNGTestRig {
 public:
  ParetoType2TestRig()
      : VectorRNGTestRig(10000, 10, {-2.5, -1.7, -0.1, 0.0, 2.0, 5.8},
                         {-3, -2, -1, 0, 2, 6}, {}, {}, {0.1, 1.0, 2.5, 4.0},
                         {1, 2, 3, 4}, {-2.7, -1.5, -0.5, 0.0}, {-3, -2, -1, 0},
                         {0.2, 1.1, 2.5, 2.0}, {1, 2, 3, 4},
                         {-3.7, -2.5, -1.0, 0.0}, {-4, -3, -2, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& mu, const T2& lambda, const T3& alpha,
                        T_rng& rng) const {
    return stan::math::pareto_type_2_rng(mu, lambda, alpha, rng);
  }

  std::vector<double> generate_quantiles(double mu, double lambda,
                                         double alpha) const {
    return std::vector<double>();
  }
};

TEST(ProbDistributionsParetoType2, errorCheck) {
  check_dist_throws_all_types(ParetoType2TestRig());
}
