#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
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
