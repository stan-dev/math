#include <test/unit/math/prim/prob/hmm_util.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <stan/math/prim/prob/hmm_latent_rng.hpp>
#include <stan/math/prim/prob/chi_square_lcdf.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(hmm_rng_test, chiSquareGoodnessFitTest) {
  // With identity transition and constant log_omegas, the sampled latent
  // states are identifcal and follow a Bernoulli distribution parameterized
  // by rho.
  // The samples live on {1, 2}, so we need to deduct the error_index to
  // to make the indices 0-indexed.
  using stan::math::hmm_latent_rng;

  int n_states = 2;
  int n_transitions = 10;
  Eigen::MatrixXd Gamma = Eigen::MatrixXd::Identity(n_states, n_states);
  Eigen::VectorXd rho(n_states);
  rho << 0.65, 0.35;
  Eigen::MatrixXd log_omegas
      = Eigen::MatrixXd::Ones(n_states, n_transitions + 1);

  boost::random::mt19937 rng;
  int N = 10000;

  std::vector<double> expected;
  expected.push_back(N * rho(0));
  expected.push_back(N * rho(1));

  std::vector<int> counts(2);
  std::vector<int> state;

  for (int i = 0; i < N; ++i) {
    state = hmm_latent_rng(log_omegas, Gamma, rho, rng);
    for (int j = 1; j < n_states; ++j)
      EXPECT_EQ(state[j], state[0]);

    ++counts[state[0] - stan::error_index::value];
  }

  assert_chi_squared(counts, expected, 1e-6);
}

TEST(hmm_rng_test, chiSquareGoodnessFitTest_symmetric) {
  // In this two states situation, the latent states are
  // symmetric, based on the observational log density,
  // and transition matrix.
  // The initial conditions introduces an asymmetry in the first
  // state. The other hidden states all have probability 0.5.
  // Note 1: the hidden states are also uncorrelated.
  // Note 2: as before, to do a chi-squared test, we subtract
  //  the error_index from hidden_state, to produce variables
  //  on {0, 1}.
  using stan::math::hmm_latent_rng;

  int n_states = 2;
  int n_transitions = 1;
  Eigen::MatrixXd Gamma(n_states, n_states);
  Gamma << 0.5, 0.5, 0.5, 0.5;
  Eigen::VectorXd rho(n_states);
  rho << 0.3, 0.7;
  Eigen::MatrixXd log_omegas
      = Eigen::MatrixXd::Ones(n_states, n_transitions + 1);

  boost::random::mt19937 rng;
  int N = 10000;

  std::vector<double> expected_0;
  expected_0.push_back(N * rho(0));
  expected_0.push_back(N * rho(1));

  std::vector<double> expected_1;
  expected_1.push_back(N * 0.5);
  expected_1.push_back(N * 0.5);

  std::vector<int> counts_0(2);
  std::vector<int> counts_1(2);
  // int product = 0;
  std::vector<int> states;
  int a = 0, c = 0;
  for (int i = 0; i < N; ++i) {
    states = hmm_latent_rng(log_omegas, Gamma, rho, rng);
    ++counts_0[states[0] - stan::error_index::value];
    ++counts_1[states[1] - stan::error_index::value];
    // product += states[0] * states[1];
    a += (states[0] == stan::error_index::value
          && states[1] == stan::error_index::value);
    c += (states[0] == 1 + stan::error_index::value
          && states[1] == stan::error_index::value);
  }

  // Test the marginal probabilities of each variable
  assert_chi_squared(counts_0, expected_0, 1e-6);
  assert_chi_squared(counts_1, expected_1, 1e-6);

  // Test for independence (0 correlation by construction).
  // By independence E(XY) = E(X)E(Y). We compute the R.H.S
  // analytically and the L.H.S numerically.
  std::vector<int> counts_xy(2);
  counts_xy[0] = a;
  counts_xy[1] = c;
  std::vector<double> expected_xy;
  expected_xy.push_back(N * rho(0) * 0.5);
  expected_xy.push_back(N * rho(1) * 0.5);
  assert_chi_squared(counts_xy, expected_xy, 1e-6);

  // DRAFT -- code for chi-squared independence test.
  // (overkill, since we have analytical prob for each cell)
  // Test that the two states are independent, using a chi squared
  // test for independence.
  // Eigen::MatrixXd Expected(n_states, (n_transitions + 1));
  // Expected << (a + b) * (a + c), (a + b) * (b + d),
  //             (c + d) * (a + c), (c + d) * (b + d);
  // Expected = Expected / N;
  //
  // Eigen::MatrixXd Observed(n_states, (n_transitions + 1));
  // Observed << a, b, c, d;
  // double chi = 0;
  //
  // for (int i = 0; i < n_states; ++i)
  //   for (int j = 0; j < n_transitions + 1; ++j)
  //     chi += (Observed(i, j) - Expected(i, j))
  //             * (Observed(i, j) - Expected(i, j)) / Expected(i, j);
  //
  // int nu = 1;
  // double p_value = exp(stan::math::chi_square_lcdf(chi, nu));
  // double threshold = 0.1;  // CHECK -- what is an appropriate threshold?
  // EXPECT_TRUE(p_value > threshold);
}
