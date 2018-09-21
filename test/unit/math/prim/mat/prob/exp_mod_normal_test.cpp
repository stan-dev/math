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
};

TEST(ProbDistributionsExpModNormal, errorCheck) {
  check_dist_throws_all_types(ExpModNormalTestRig());
}

/*
 * Don't have an easy way to compute ExpModNormal quantiles in C++, so test
 * the distributions manually
 */
TEST(ProbDistributionsExpModNormal, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;

  double mu = 2.0;
  double sigma = 1.0;
  double lambda = 3.0;

  std::vector<double> loc
      = {0.00508630279050534, 0.2937255454905956, 0.478096389341816,
         0.6175435407098293,  0.7315151474840027, 0.8289498399078623,
         0.9147360619484722,  0.9918540497848275, 1.062261803414033,
         1.127318765520224,   1.188009350988437,  1.245070983987459,
         1.299071267620822,   1.350457156638078,  1.399587697999158,
         1.446755879713708,   1.492204484340773,  1.536137286203502,
         1.578727150279017,   1.620122502460362,  1.660451793116215,
         1.699827269095549,   1.738347835774777,  1.776101227348508,
         1.813165978184319,   1.849612768859248,  1.885505736729774,
         1.920903407354275,   1.955859553768006,  1.990423920241132,
         2.024642826361759,   2.05855965364628,   2.092215318521506,
         2.125648663501404,   2.158896791954166,  2.191995413723774,
         2.224979095736312,   2.257881543752706,  2.290735857678737,
         2.323574747626378,   2.356430775482149,  2.389336578988231,
         2.422325095366644,   2.455429795617495,  2.488684916286868,
         2.522125743327702,   2.555788838498307,  2.589712378661565,
         2.623936478994747,   2.658503537633071,  2.693458675183349,
         2.728850235921933,   2.764730264946303,  2.801155243664606,
         2.838186780108561,   2.875892548595454,  2.91434736079305,
         2.953634467630173,   2.993847166352555,  3.035090766145558,
         3.077485080687957,   3.121167522841614,  3.166296968905866,
         3.213059109917944,   3.261673003998312,  3.312400228962045,
         3.365557130017247,   3.421531936000286,  3.480809161629723,
         3.544005377092386,   3.611923435318303,  3.685638028101141,
         3.766637605625211,   3.857072937613304,  3.960230401442832,
         4.081526033985821,   4.230921381548062,  4.430265245788379,
         4.747088207522281,   10.17028381249956};
  int K = loc.size();
  boost::math::chi_squared mydist(K - 1);

  std::vector<double> a1 = stan::math::exp_mod_normal_rng(
      mu, std::vector<double>(N, sigma), lambda, rng);
  assert_matches_quantiles(a1, loc, 1e-6);
  std::vector<double> a2(N, 0.0);
  for (size_t i = 0; i < a2.size(); ++i) {
    std::vector<double> mu_ = {7.0, mu};
    std::vector<double> sigma_ = {2.0, sigma};
    std::vector<double> lambda_ = {1.0, lambda};
    a2[i] = stan::math::exp_mod_normal_rng(mu_, sigma_, lambda_, rng)[1];
  }
  assert_matches_quantiles(a2, loc, 1e-6);
  std::vector<double> a3(N, 0.0);
  for (size_t i = 0; i < a2.size(); ++i) {
    Eigen::VectorXd mu_(2);
    mu_ << 7.0, mu;
    Eigen::RowVectorXd sigma_(2);
    sigma_ << 1.0, sigma;
    Eigen::RowVectorXd lambda_(2);
    lambda_ << 1.0, lambda;
    a3[i] = stan::math::exp_mod_normal_rng(mu_, sigma_, lambda_, rng)[1];
  }
  assert_matches_quantiles(a3, loc, 1e-6);
}
