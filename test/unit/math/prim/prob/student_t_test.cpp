#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

class StudentTTestRig : public VectorRealRNGTestRig {
 public:
  StudentTTestRig()
      : VectorRealRNGTestRig(10000, 10, {1.1, 2.0, 2.5, 3.1}, {1, 2, 3, 4},
                             {-1.7, -0.5, -2.5, 0.0}, {-2, -1, -3, 0},
                             {-2.5, -1.7, -0.1, 0.0, 2.0, 5.8},
                             {-3, -2, -1, 0, 2, 6}, {}, {},
                             {0.1, 1.0, 2.5, 4.0}, {1, 2, 3, 4},
                             {-2.7, -1.5, -0.5, 0.0}, {-3, -2, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& nu, const T2& mu, const T3& sigma,
                        T_rng& rng) const {
    return stan::math::student_t_rng(nu, mu, sigma, rng);
  }

  std::vector<double> generate_quantiles(double nu, double mu,
                                         double sigma) const {
    std::vector<double> quantiles;
    double K = stan::math::round(2 * std::pow(N_, 0.4));

    boost::math::students_t_distribution<> dist(nu);
    for (int i = 1; i < K; ++i) {
      double frac = i / K;
      quantiles.push_back(quantile(dist, frac) * sigma + mu);
    }
    quantiles.push_back(std::numeric_limits<double>::max());

    return quantiles;
  }
};

TEST(ProbDistributionsStudentT, errorCheck) {
  check_dist_throws_all_types(StudentTTestRig());
}

TEST(ProbDistributionsStudentT, distributionTest) {
  check_quantiles_real_real_real(StudentTTestRig());
}

TEST(ProbDistributionsStudentT, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::student_t_rng(3.0, 2.0, 2.0, rng));

  EXPECT_THROW(stan::math::student_t_rng(3.0, 2.0, -2.0, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::student_t_rng(-3.0, 2.0, 2.0, rng),
               std::domain_error);
  EXPECT_THROW(
      stan::math::student_t_rng(stan::math::positive_infinity(), 2.0, 2.0, rng),
      std::domain_error);
  EXPECT_THROW(
      stan::math::student_t_rng(3, stan::math::positive_infinity(), 2.0, rng),
      std::domain_error);
  EXPECT_THROW(
      stan::math::student_t_rng(3, 2, stan::math::positive_infinity(), rng),
      std::domain_error);
}

TEST(ProbDistributionsStudentT, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));

  std::vector<double> samples;
  for (int i = 0; i < N; ++i) {
    samples.push_back((stan::math::student_t_rng(3.0, 2.0, 2.0, rng) - 2.0)
                      / 2.0);
  }

  // Generate quantiles from boost's student t distribution
  boost::math::students_t_distribution<> dist(3.0);
  std::vector<double> quantiles;
  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(quantile(dist, frac));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  // Assert that they match
  assert_matches_quantiles(samples, quantiles, 1e-6);
}
