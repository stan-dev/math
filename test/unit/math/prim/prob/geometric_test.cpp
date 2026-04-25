#include <stan/math/prim.hpp>
#include <stan/math/prim/prob/geometric_rng.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/VectorIntRNGTestRig.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <limits>
#include <vector>

class GeometricTestRig : public VectorIntRNGTestRig {
 public:
  GeometricTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                            {0.1, 0.3, 0.5, 0.7, 0.9}, {1}, {-0.1, 1.1},
                            {-1, 0, 2}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& theta, const T2&, const T3&,
                        T_rng& rng) const {
    return stan::math::geometric_rng(theta, rng);
  }

  template <typename T1>
  double pmf(int y, T1 theta, double, double) const {
    if (y < 0) {
      return 0.0;
    }
    return std::exp(stan::math::geometric_lpmf(y, theta));
  }
};

TEST(ProbDistributionsGeometric, errorCheck) {
  check_dist_throws_all_types(GeometricTestRig());
}

TEST(ProbDistributionsGeometric, distributionCheck) {
  check_counts_real(GeometricTestRig());
}

TEST(ProbDistributionsGeometric, error_check) {
  boost::random::mt19937 rng;

  EXPECT_NO_THROW(stan::math::geometric_rng(0.5, rng));
  EXPECT_NO_THROW(stan::math::geometric_rng(0.1, rng));
  EXPECT_NO_THROW(stan::math::geometric_rng(0.9, rng));
  EXPECT_NO_THROW(stan::math::geometric_rng(1.0, rng));

  EXPECT_THROW(stan::math::geometric_rng(0.0, rng), std::domain_error);
  EXPECT_THROW(stan::math::geometric_rng(-0.5, rng), std::domain_error);

  EXPECT_THROW(stan::math::geometric_rng(stan::math::positive_infinity(), rng),
               std::domain_error);
  EXPECT_THROW(stan::math::geometric_rng(stan::math::not_a_number(), rng),
               std::domain_error);
}

TEST(ProbDistributionsGeometric, lpmf_values) {
  // P(n=0|theta=0.5) = 0.5
  EXPECT_NEAR(stan::math::geometric_lpmf(0, 0.5), std::log(0.5), 1e-10);
  // P(n=2|theta=0.5) = 0.5 * 0.5^2 = 0.125
  EXPECT_NEAR(stan::math::geometric_lpmf(2, 0.5), std::log(0.125), 1e-10);
  // P(n=0|theta=1.0) = 1.0
  EXPECT_NEAR(stan::math::geometric_lpmf(0, 1.0), 0.0, 1e-10);
  // P(n=1|theta=1.0) = 0
  EXPECT_EQ(stan::math::geometric_lpmf(1, 1.0),
            stan::math::negative_infinity());
  // P(n=0|theta=0.3) = 0.3
  EXPECT_NEAR(stan::math::geometric_lpmf(0, 0.3), std::log(0.3), 1e-10);
  // P(n=3|theta=0.3) = 0.3 * 0.7^3 = 0.1029
  EXPECT_NEAR(stan::math::geometric_lpmf(3, 0.3), std::log(0.1029), 1e-10);
}

TEST(ProbDistributionsGeometric, cdf_values) {
  // CDF(n=0|theta=0.5) = 1 - (1-0.5)^1 = 0.5
  EXPECT_NEAR(stan::math::geometric_cdf(0, 0.5), 0.5, 1e-10);
  // CDF(n=2|theta=0.5) = 1 - (0.5)^3 = 0.875
  EXPECT_NEAR(stan::math::geometric_cdf(2, 0.5), 0.875, 1e-10);
  // CDF(n=0|theta=1.0) = 1.0
  EXPECT_NEAR(stan::math::geometric_cdf(0, 1.0), 1.0, 1e-10);
  // CDF(n=4|theta=0.3) = 1 - 0.7^5 = 0.83193
  EXPECT_NEAR(stan::math::geometric_cdf(4, 0.3), 0.83193, 1e-5);
}

TEST(ProbDistributionsGeometric, cdf_below_support) {
  EXPECT_DOUBLE_EQ(stan::math::geometric_cdf(-1, 0.5), 0.0);
  EXPECT_DOUBLE_EQ(stan::math::geometric_lcdf(-1, 0.5),
                   stan::math::negative_infinity());
}

TEST(ProbDistributionsGeometric, lccdf_values) {
  // CCDF(n=0|theta=0.5) = (1-0.5)^1 = 0.5
  EXPECT_NEAR(stan::math::geometric_lccdf(0, 0.5), std::log(0.5), 1e-10);
  // CCDF(n=2|theta=0.5) = (0.5)^3 = 0.125
  EXPECT_NEAR(stan::math::geometric_lccdf(2, 0.5), std::log(0.125), 1e-10);
}

TEST(ProbDistributionsGeometric, rng_theta_one) {
  boost::random::mt19937 rng;
  // theta=1.0 should always return 0 (always succeed on first trial)
  for (int i = 0; i < 100; i++) {
    EXPECT_EQ(stan::math::geometric_rng(1.0, rng), 0);
  }
}

TEST(ProbDistributionsGeometric, theta_zero_throws) {
  // theta = 0 is a degenerate parameter (no successes possible) and is
  // outside the documented domain (0, 1]; all four distribution functions
  // must throw rather than returning a poisoned partial.
  EXPECT_THROW(stan::math::geometric_lpmf(0, 0.0), std::domain_error);
  EXPECT_THROW(stan::math::geometric_lpmf(3, 0.0), std::domain_error);
  EXPECT_THROW(stan::math::geometric_cdf(0, 0.0), std::domain_error);
  EXPECT_THROW(stan::math::geometric_lcdf(0, 0.0), std::domain_error);
  EXPECT_THROW(stan::math::geometric_lccdf(0, 0.0), std::domain_error);
}
