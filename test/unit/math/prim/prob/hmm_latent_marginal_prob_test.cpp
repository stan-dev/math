#include <test/unit/math/prim/prob/hmm_marginal_test.cpp>
#include <stan/math/prim/prob/hmm_latent_marginal_prob.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>


TEST_F(hmm_test, latent_prob_transition_1) {
  using stan::math::hmm_latent_marginal_prob;

  int n_states = 2;
  Eigen::MatrixXd Gamma(n_states, n_states);
  Gamma << 1, 0, 1, 0;
  Eigen::VectorXd rho(n_states);
  rho << 1, 0;

  Eigen::MatrixXd prob = hmm_latent_marginal_prob(log_omegas_, Gamma, rho);

  for (int i = 0; i < n_transitions_; i++) EXPECT_EQ(prob(0, i), 1);
  for (int i = 0; i < n_transitions_; i++) EXPECT_EQ(prob(1, i), 0);
}
