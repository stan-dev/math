#include <test/unit/math/prim/prob/hmm_marginal_test.cpp>
#include <stan/math/prim/prob/hmm_latent_marginal_prob.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST_F(hmm_test, latent_prob_single_outcome) {
  using stan::math::hmm_latent_marginal_prob;

  int n_states = 2;
  Eigen::MatrixXd Gamma(n_states, n_states);
  Gamma << 1, 0, 1, 0;
  Eigen::VectorXd rho(n_states);
  rho << 1, 0;

  Eigen::MatrixXd prob = hmm_latent_marginal_prob(log_omegas_, Gamma, rho);

  for (int i = 0; i < n_transitions_; i++) {
    EXPECT_EQ(prob(0, i), 1);
    EXPECT_EQ(prob(1, i), 0);
  }
}

TEST_F(hmm_test, latent_prob_identity_transition) {
  // With an identity transition matrix, all latent probabilities
  // are equal. Setting the log density to 1 for all states makes
  // the initial prob drive the subsequent probabilities.
  using stan::math::hmm_latent_marginal_prob;
  int n_states = 2;
  Eigen::MatrixXd Gamma = Eigen::MatrixXd::Identity(n_states, n_states);
  Eigen::MatrixXd log_omegas
      = Eigen::MatrixXd::Ones(n_states, n_transitions_ + 1);

  Eigen::MatrixXd prob = hmm_latent_marginal_prob(log_omegas, Gamma, rho_);

  for (int i = 0; i < n_transitions_; i++) {
    EXPECT_FLOAT_EQ(prob(0, i), rho_(0));
    EXPECT_FLOAT_EQ(prob(1, i), rho_(1));
  }
}
