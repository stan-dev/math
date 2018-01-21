#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <vector>

class ExpModNormalTestRig : public VectorRNGTestRig {
 public:
  ExpModNormalTestRig()
      : VectorRNGTestRig(10000, 10, {-2.5, -1.7, -0.1, 0.0, 2.0, 5.8},
                         {-3, -2, -1, 0, 2, 6}, {}, {}, {0.1, 1.0, 2.5, 4.0},
                         {1, 2, 3, 4}, {-2.7, -1.5, -0.5, 0.0}, {-3, -2, -1, 0},
                         {0.2, 1.1, 2.5, 2.0}, {1, 2, 3, 4},
                         {-3.7, -2.5, -1.0, 0.0}, {-4, -3, -2, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& loc, const T2& scale, const T3& inv_scale,
                        T_rng& rng) const {
    return stan::math::exp_mod_normal_rng(loc, scale, inv_scale, rng);
  }

  std::vector<double> generate_quantiles(double loc, double scale,
                                         double inv_scale) const {
    return std::vector<double>();
  }
};

TEST(ProbDistributionsExpModNormal, errorCheck) {
  check_dist_throws_all_types(ExpModNormalTestRig());
}
