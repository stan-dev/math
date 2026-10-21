#include <test/unit/math/prim/prob/hmm_util.hpp>
#include <stan/math/prim/prob/hmm_hidden_state_prob.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

namespace {
/**
 * Reference implementation of hmm_hidden_state_prob: the forward-backward
 * algorithm computed entirely on the log scale, with log_sum_exp in place
 * of every sum, so that neither the forward nor the backward pass can
 * underflow or overflow.
 */
Eigen::MatrixXd hmm_hidden_state_prob_log_space(
    const Eigen::MatrixXd& log_omegas, const Eigen::MatrixXd& Gamma,
    const Eigen::VectorXd& rho) {
  using stan::math::log_sum_exp;
  const int n_states = log_omegas.rows();
  const int n_cols = log_omegas.cols();
  const Eigen::MatrixXd log_Gamma = Gamma.array().log();

  Eigen::MatrixXd log_alphas(n_states, n_cols);
  log_alphas.col(0) = rho.array().log() + log_omegas.col(0).array();
  for (int n = 1; n < n_cols; ++n) {
    for (int j = 0; j < n_states; ++j) {
      Eigen::VectorXd terms = log_alphas.col(n - 1) + log_Gamma.col(j);
      log_alphas(j, n) = log_omegas(j, n) + log_sum_exp(terms);
    }
  }

  Eigen::MatrixXd prob(n_states, n_cols);
  Eigen::VectorXd log_beta = Eigen::VectorXd::Zero(n_states);
  for (int n = n_cols - 1; n >= 0; --n) {
    if (n < n_cols - 1) {
      Eigen::VectorXd next = log_omegas.col(n + 1) + log_beta;
      for (int i = 0; i < n_states; ++i) {
        Eigen::VectorXd terms = log_Gamma.row(i).transpose() + next;
        log_beta(i) = log_sum_exp(terms);
      }
    }
    Eigen::VectorXd log_prob = log_alphas.col(n) + log_beta;
    prob.col(n) = (log_prob.array() - log_sum_exp(log_prob)).exp();
  }
  return prob;
}
}  // namespace

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

TEST(hmm_test_nonstandard, hidden_state_prob1) {
  // This time, the transition matrix forces states to transition
  // to the first state with probability 1. The large log densities
  // overflow if they are exponentiated directly (#2677).
  using stan::math::hmm_hidden_state_prob;
  int n_states = 2;
  int n_transitions = 2;
  Eigen::MatrixXd Gamma(n_states, n_states);
  Gamma << 1, 0, 1, 0;
  Eigen::VectorXd rho(n_states);
  rho << 0.3, 0.7;
  Eigen::MatrixXd log_omegas
      = 1000 * Eigen::MatrixXd::Ones(n_states, n_transitions + 1);

  Eigen::MatrixXd prob = hmm_hidden_state_prob(log_omegas, Gamma, rho);
  EXPECT_FLOAT_EQ(prob(0, 0), 0.3);
  EXPECT_FLOAT_EQ(prob(1, 0), 0.7);

  for (int i = 1; i < n_transitions; i++) {
    EXPECT_FLOAT_EQ(prob(0, i), 1);
    EXPECT_FLOAT_EQ(prob(1, i), 0);
  }
}

TEST_F(hmm_test, hidden_state_very_negative_log_omegas) {
  // Regression test for #2677. Every log density is below -1000, so
  // exponentiating it directly underflows to 0. Adding a constant to a
  // column of log_omegas does not change the hidden state probabilities.
  using stan::math::hmm_hidden_state_prob;
  Eigen::MatrixXd log_omegas = log_omegas_;
  for (int n = 0; n < log_omegas.cols(); ++n) {
    log_omegas.col(n).array() -= 1000 + 10 * n;
  }

  Eigen::MatrixXd prob = hmm_hidden_state_prob(log_omegas, Gamma_, rho_);
  Eigen::MatrixXd prob_unshifted
      = hmm_hidden_state_prob(log_omegas_, Gamma_, rho_);
  Eigen::MatrixXd prob_ref
      = hmm_hidden_state_prob_log_space(log_omegas, Gamma_, rho_);

  for (int n = 0; n < prob.cols(); ++n) {
    for (int k = 0; k < prob.rows(); ++k) {
      EXPECT_TRUE(std::isfinite(prob(k, n)));
      EXPECT_NEAR(prob(k, n), prob_unshifted(k, n), 1e-12);
      // The log scale reference loses about machine epsilon times the
      // magnitude of the log densities, hence the looser tolerance.
      EXPECT_NEAR(prob(k, n), prob_ref(k, n), 1e-10);
    }
  }
}

TEST(hmm_test_nonstandard, hidden_state_long_sequence) {
  // Regression test for #2677. The log density of a long sequence is far
  // below the smallest double, so the forward pass underflows to 0 unless
  // it is rescaled at every step.
  using stan::math::hmm_hidden_state_prob;
  using stan::math::hmm_marginal;
  int n_states = 2;
  int n_obs = 5000;
  Eigen::MatrixXd Gamma(n_states, n_states);
  Gamma << 0.95, 0.05, 0.10, 0.90;
  Eigen::VectorXd rho(n_states);
  rho << 2.0 / 3, 1.0 / 3;

  // Simulate a two-state Gaussian hidden Markov model, with observations
  // normal(1, 1) in state 0 and normal(-1, 1) in state 1.
  boost::random::mt19937 rng(2677);
  boost::random::uniform_01<> uniform;
  boost::random::normal_distribution<> normal(0, 1);
  Eigen::MatrixXd log_omegas(n_states, n_obs);
  int state = 0;
  for (int n = 0; n < n_obs; ++n) {
    state = uniform(rng) < Gamma(state, 0) ? 0 : 1;
    double y = (state == 0 ? 1 : -1) + normal(rng);
    for (int k = 0; k < n_states; ++k) {
      log_omegas(k, n) = state_lpdf(y, 1, 1, k);
    }
  }
  // The joint density underflows the double range many times over.
  EXPECT_LT(hmm_marginal(log_omegas, Gamma, rho), -2000);

  Eigen::MatrixXd prob = hmm_hidden_state_prob(log_omegas, Gamma, rho);
  Eigen::MatrixXd prob_ref
      = hmm_hidden_state_prob_log_space(log_omegas, Gamma, rho);

  int first_non_finite = -1;
  for (int n = 0; n < n_obs && first_non_finite < 0; ++n) {
    if (!prob.col(n).allFinite()) {
      first_non_finite = n;
    }
  }
  ASSERT_EQ(first_non_finite, -1) << "first non-finite column";
  EXPECT_LT((prob.colwise().sum().array() - 1).abs().maxCoeff(), 1e-12);
  // The log scale reference loses about machine epsilon times the
  // magnitude of the log densities, hence the looser tolerance.
  EXPECT_LT((prob - prob_ref).cwiseAbs().maxCoeff(), 1e-10);
}
