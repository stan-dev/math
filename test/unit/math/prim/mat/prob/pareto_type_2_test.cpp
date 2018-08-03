#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/prob/vector_rng_test_helper.hpp>
#include <vector>
#include <limits>

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
};

TEST(ProbDistributionsParetoType2, errorCheck) {
  check_dist_throws_all_types(ParetoType2TestRig());
}

/*
 * Don't have an easy way to compute Pareto Type 2 quantiles in C++, so test
 * the distributions manually
 */
TEST(ProbDistributionsParetoType2, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;

  double mu = 3.0;
  double lambda = 2.0;
  double alpha = 4.0;

  std::vector<double> loc
      = {3.006299290587777672812, 3.012699050612054563203,
         3.019202201737360535816, 3.025811789959920172777,
         3.032530992618458931531, 3.039363125899919637618,
         3.046311652882255582853, 3.053380192160681794178,
         3.060572527108473828150, 3.067892615828682334467,
         3.075344601859007820366, 3.082932825698688539262,
         3.090661837233660502733, 3.098536409144583281261,
         3.106561551391706199610, 3.114742526881128181770,
         3.123084868428977589616, 3.131594397153576458948,
         3.140277242441022131914, 3.149139863647083981135,
         3.158189073718207673380, 3.167432064937129254645,
         3.176876437024558796907, 3.186530227858186847811,
         3.196401947104449803305, 3.206500613097915231009,
         3.216835793348621130860, 3.227417649110363662146,
         3.238256984504057189866, 3.249365300761396291307,
         3.260754856237075394176, 3.272438732934998650848,
         3.284430910408003168044, 3.296746348025026307482,
         3.309401076758502924235, 3.322412301833146841545,
         3.335798517801522145021, 3.349579637879685378721,
         3.363777139697528006934, 3.378414230005442053795,
         3.393516031347077355207, 3.409109794276735438956,
         3.425225139394374007651, 3.441894334323137627507,
         3.459152611805057286887, 3.477038536393603163077,
         3.495594428846999779381, 3.514866859365870954690,
         3.534907223394969211228, 3.555772416985089723340,
         3.577525632908273856003, 3.600237304137476890276,
         3.623986228353907268485, 3.648860916404678800262,
         3.674961219905688025733, 3.702400309614068518727,
         3.731307098427877555480, 3.761829233309825326614,
         3.794136822632038885672, 3.828427124746190290949,
         3.864930508176635193252, 3.903918116461908560666,
         3.945711851615199439181, 3.990697562442441448383,
         4.039342742606370606495, 4.092220702775781759897,
         4.150044248441071559341, 4.213713675778606670974,
         4.284386974674554871001, 4.363585661014857564055,
         4.453360085481801178275, 4.556558820077846405638,
         4.677290860899674385109, 4.821771168817467589918,
         5.000000000000000000000, 5.229485053762255475363,
         5.544877465869998189874, 6.029733718731740133023,
         6.981395124884887337657, std::numeric_limits<double>::infinity()};

  int K = loc.size();
  boost::math::chi_squared mydist(K - 1);

  std::vector<double> a1 = stan::math::pareto_type_2_rng(
      mu, lambda, std::vector<double>(N, alpha), rng);
  assert_matches_quantiles(a1, loc, 1e-6);
  std::vector<double> a2(N, 0.0);
  for (size_t i = 0; i < a2.size(); ++i) {
    std::vector<double> mu_ = {7.0, mu};
    std::vector<double> lambda_ = {2.0, lambda};
    std::vector<double> alpha_ = {1.0, alpha};
    a2[i] = stan::math::pareto_type_2_rng(mu_, lambda_, alpha_, rng)[1];
  }
  assert_matches_quantiles(a2, loc, 1e-6);
  std::vector<double> a3(N, 0.0);
  for (size_t i = 0; i < a2.size(); ++i) {
    Eigen::VectorXd mu_(2);
    mu_ << 7.0, mu;
    Eigen::RowVectorXd lambda_(2);
    lambda_ << 1.0, lambda;
    Eigen::RowVectorXd alpha_(2);
    alpha_ << 1.0, alpha;
    a3[i] = stan::math::pareto_type_2_rng(mu_, lambda_, alpha_, rng)[1];
  }
  assert_matches_quantiles(a3, loc, 1e-6);
}
