#include <stan/math/prim.hpp>
#include <stan/math/prim/prob/yule_simon_rng.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/VectorIntRNGTestRig.hpp>
#include <gtest/gtest.h>
#include <random>
#include <boost/math/distributions.hpp>
#include <limits>
#include <vector>

class YuleSimonTestRig : public VectorIntRNGTestRig {
 public:
  YuleSimonTestRig()
      : VectorIntRNGTestRig(10000, 10, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                            {2.1, 4.1, 10.1}, {1, 2, 3}, {-3.0, -2.0, 0.0},
                            {-3, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& alpha, const T2&, const T3&,
                        T_rng& rng) const {
    return stan::math::yule_simon_rng(alpha, rng);
  }

  template <typename T1>
  double pmf(int y, T1 alpha, double, double) const {
    return std::exp(stan::math::yule_simon_lpmf(y, alpha));
  }
};

TEST(ProbDistributionsYuleSimon, errorCheck) {
  check_dist_throws_all_types(YuleSimonTestRig());
}

TEST(ProbDistributionsYuleSimon, distributionCheck) {
  check_counts_real(YuleSimonTestRig());
}

TEST(ProbDistributionsYuleSimon, error_check) {
  std::mt19937 rng;

  EXPECT_NO_THROW(stan::math::yule_simon_rng(1.0, rng));
  EXPECT_NO_THROW(stan::math::yule_simon_rng(2.0, rng));
  EXPECT_NO_THROW(stan::math::yule_simon_rng(10.0, rng));

  EXPECT_THROW(stan::math::yule_simon_rng(-0.5, rng), std::domain_error);
  EXPECT_THROW(stan::math::yule_simon_rng(0.0, rng), std::domain_error);
  EXPECT_THROW(stan::math::yule_simon_rng(-10.0, rng), std::domain_error);

  EXPECT_THROW(stan::math::yule_simon_rng(stan::math::positive_infinity(), rng),
               std::domain_error);

  EXPECT_THROW(stan::math::yule_simon_rng(stan::math::not_a_number(), rng),
               std::domain_error);
}
