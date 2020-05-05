#include <test/unit/math/prim/prob/hmm_marginal_test.cpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <stan/math/prim/prob/hmm_latent_rng.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>


// TEST_F(hmm_test, compile_and_run) {
//   using stan::math::hmm_latent_rng;
//   boost::random::mt19937 rng;
//
//   std::vector<int> x = hmm_latent_rng(log_omegas_, Gamma_, rho_, rng);
//   for (size_t i = 0; i < x.size(); ++i)
//     std::cout << x[i] << " ";
//   std::cout << std::endl;
// }

TEST(hmm_rng_test, chiSquareGoodnessFitTest) {
  // with identity transition and constant log_omegas, the sampled latent
  // states are identifcal and follow a Bernoulli distribution parameterized
  // by rho.
  using stan::math::hmm_latent_rng;

  int n_states = 2;
  int n_transitions = 10;
  Eigen::MatrixXd Gamma = Eigen::MatrixXd::Identity(n_states, n_states);
  Eigen::VectorXd rho(n_states);
  rho << 0.65, 0.35;
  Eigen::MatrixXd log_omegas =
    Eigen::MatrixXd::Ones(n_states, n_transitions + 1);

  boost::random::mt19937 rng;
  int N = 10000;

  std::vector<double> expected;
  expected.push_back(N * rho(0));
  expected.push_back(N * rho(1));

  std::vector<int> counts(2);
  std::vector<int> state;
  for (int i = 0; i < N; ++i) {
    state = hmm_latent_rng(log_omegas, Gamma, rho, rng);
    for (int j = 1; j < n_states; ++j) EXPECT_EQ(state[j], state[0]);
    ++counts[state[0]];
  }

  assert_chi_squared(counts, expected, 1e-6);
}
