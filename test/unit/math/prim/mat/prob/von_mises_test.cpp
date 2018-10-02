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
      : VectorRNGTestRig(10000, 10, {-2.5, -1.7, -0.1, 0.0, 2.0, 5.8},
                         {-3, -2, -1, 0, 2, 6}, {}, {}, {0.1, 1.0, 2.5, 4.0},
                         {1, 2, 3, 4}, {-2.7, -1.5, -0.5, 0.0},
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

/*
 * Don't have an easy way to compute VonMises quantiles in C++, so test
 * the distributions manually
 */
TEST(ProbDistributionsVonMises, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;

  std::vector<double> loc
      = {1.589032, 1.835082,         1.977091, 2.078807, 2.158994, 2.225759,
         2.283346, 2.334261,         2.380107, 2.421973, 2.460635, 2.496664,
         2.530492, 2.562457,         2.592827, 2.621819, 2.649609, 2.676346,
         2.702153, 2.727137,         2.751389, 2.774988, 2.798002, 2.820493,
         2.842515, 2.864116,         2.885340, 2.906227, 2.926814, 2.947132,
         2.967214, 2.987089,         3.006782, 3.026321, 3.045728, 3.065028,
         3.084243, 3.103394,         3.122504, 3.141593, 3.160681, 3.179791,
         3.198942, 3.218157,         3.237457, 3.256865, 3.276403, 3.296097,
         3.315971, 3.336053,         3.356372, 3.376958, 3.397845, 3.419069,
         3.440671, 3.462693,         3.485184, 3.508198, 3.531796, 3.556048,
         3.581032, 3.606840,         3.633576, 3.661367, 3.690358, 3.720728,
         3.752694, 3.786522,         3.822550, 3.861212, 3.903079, 3.948925,
         3.999839, 4.057427,         4.124191, 4.204379, 4.306094, 4.448103,
         4.694153, 6.283185307179586};
  int K = loc.size();
  boost::math::chi_squared mydist(K - 1);
  for (int i = 0; i < K; i++)
    loc[i] = loc[i] - stan::math::pi();

  std::vector<double> a1
      = stan::math::von_mises_rng(0, std::vector<double>(N, 3.0), rng);
  assert_matches_quantiles(a1, loc, 1e-6);
  std::vector<double> a2(N, 0.0);
  for (size_t i = 0; i < a2.size(); ++i) {
    std::vector<double> mu = {7.0, 0.0};
    std::vector<double> kappa = {1.0, 3.0};
    a2[i] = stan::math::von_mises_rng(mu, kappa, rng)[1];
  }
  assert_matches_quantiles(a2, loc, 1e-6);
  std::vector<double> a3(N, 0.0);
  for (size_t i = 0; i < a2.size(); ++i) {
    Eigen::VectorXd mu(2);
    mu << 7.0, 0.0;
    Eigen::RowVectorXd kappa(2);
    kappa << 1.0, 3.0;
    a3[i] = stan::math::von_mises_rng(mu, kappa, rng)[1];
  }
  assert_matches_quantiles(a3, loc, 1e-6);
}
