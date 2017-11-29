#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/mat/prob/VectorRNGTestRig.hpp>
#include <limits>
#include <vector>

class ChiSquareTestRig : public VectorRNGTestRig {
public:
  ChiSquareTestRig() :
    VectorRNGTestRig(10000, 10,
                     {0.1, 1.0, 2.5, 4.0},
                     {1, 2, 3, 4},
                     {-2.7, -1.5, -0.5, 0.0},
                     {-3, -2, -1, 0}) {}

  template<typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& nu, const T2&, const T3&,
                        T_rng& rng) const {
    return stan::math::chi_square_rng(nu, rng);
  }

  std::vector<double> generate_quantiles(double nu, double, double)
    const {
    std::vector<double> quantiles;
    double K = boost::math::round(2 * std::pow(N_, 0.4));
    boost::math::chi_squared_distribution<> dist(nu);

    for (int i = 1; i < K; ++i) {
      double frac = i / K;
      quantiles.push_back(quantile(dist, frac));
    }
    quantiles.push_back(std::numeric_limits<double>::max());

    return quantiles;
  }
};

TEST(ProbDistributionsChiSquare, errorCheck) {
  check_dist_throws_all_types(ChiSquareTestRig());
}

TEST(ProbDistributionsChiSquare, chiSquareGoodnessFitTest) {
  check_quantiles_all_types(ChiSquareTestRig());
}
