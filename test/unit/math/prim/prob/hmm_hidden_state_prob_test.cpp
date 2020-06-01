#include <test/unit/math/prim/prob/hmm_util.hpp>
#include <stan/math/prim/prob/hmm_hidden_state_prob.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST_F(hmm_test, hidden_state_single_outcome) {
  using stan::math::hmm_hidden_state_prob;

  int n_states = 2;
  Eigen::MatrixXd Gamma(n_states, n_states);
  Gamma << 1, 0, 1, 0;
  Eigen::VectorXd rho(n_states);
  rho << 1, 0;

  Eigen::MatrixXd prob = hmm_hidden_state_prob(log_omegas_, Gamma, rho);

  for (int i = 0; i < n_transitions_; i++) {
    EXPECT_EQ(prob(0, i), 1);
    EXPECT_EQ(prob(1, i), 0);
  }
}

TEST_F(hmm_test, hidden_state_identity_transition) {
  // With an identity transition matrix, all latent probabilities
  // are equal. Setting the log density to 1 for all states makes
  // the initial prob drive the subsequent probabilities.
  using stan::math::hmm_hidden_state_prob;
  int n_states = 2;
  Eigen::MatrixXd Gamma = Eigen::MatrixXd::Identity(n_states, n_states);
  Eigen::MatrixXd log_omegas
      = Eigen::MatrixXd::Ones(n_states, n_transitions_ + 1);

  Eigen::MatrixXd prob = hmm_hidden_state_prob(log_omegas, Gamma, rho_);

  for (int i = 0; i < n_transitions_; i++) {
    EXPECT_FLOAT_EQ(prob(0, i), rho_(0));
    EXPECT_FLOAT_EQ(prob(1, i), rho_(1));
  }
}

TEST(hmm_test_nonstandard, hidden_state_symmetry) {
  // In this two states situation, the latent states are
  // symmetric, based on the observational log density,
  // and transition matrix.
  // The initial conditions introduces an asymmetry in the first
  // state. The other hidden states all have probability 0.5.
  using stan::math::hmm_hidden_state_prob;
  int n_states = 2;
  int n_transitions = 2;
  Eigen::MatrixXd Gamma(n_states, n_states);
  Gamma << 0.5, 0.5, 0.5, 0.5;
  Eigen::VectorXd rho(n_states);
  rho << 0.3, 0.7;
  Eigen::MatrixXd log_omegas
      = Eigen::MatrixXd::Ones(n_states, n_transitions + 1);

  Eigen::MatrixXd prob = hmm_hidden_state_prob(log_omegas, Gamma, rho);

  EXPECT_FLOAT_EQ(prob(0, 0), 0.3);
  EXPECT_FLOAT_EQ(prob(1, 0), 0.7);

  for (int i = 1; i < n_transitions; i++) {
    EXPECT_FLOAT_EQ(prob(0, i), 0.5);
    EXPECT_FLOAT_EQ(prob(1, i), 0.5);
  }
}
