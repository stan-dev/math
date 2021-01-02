#ifndef TEST_UNIT_MATH_PRIM_PROB_HMM_UTIL
#define TEST_UNIT_MATH_PRIM_PROB_HMM_UTIL
#include <stan/math/prim/prob/hmm_marginal.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

/**
 * Wrapper around hmm_marginal_density which passes rho and
 * Gamma without the last element of each column. We recover
 * the last element using the fact each column sums to 1.
 * The purpose of this function is to do finite diff benchmarking,
 * without breaking the simplex constraint.
 */
template <typename T_omega, typename T_Gamma, typename T_rho>
inline stan::return_type_t<T_omega, T_Gamma, T_rho> hmm_marginal_test_wrapper(
    const Eigen::Matrix<T_omega, Eigen::Dynamic, Eigen::Dynamic>& log_omegas,
    const Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic>&
        Gamma_unconstrained,
    const std::vector<T_rho>& rho_unconstrained) {
  using stan::math::row;
  using stan::math::sum;
  int n_states = log_omegas.rows();

  Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic> Gamma(n_states,
                                                               n_states);
  for (int i = 0; i < n_states; i++) {
    Gamma(i, n_states - 1) = 1 - sum(row(Gamma_unconstrained, i + 1));
    for (int j = 0; j < n_states - 1; j++) {
      Gamma(i, j) = Gamma_unconstrained(i, j);
    }
  }

  Eigen::Matrix<T_rho, Eigen::Dynamic, 1> rho(n_states);
  rho(1) = 1 - sum(rho_unconstrained);
  for (int i = 0; i < n_states - 1; i++)
    rho(i) = rho_unconstrained[i];

  return stan::math::hmm_marginal(log_omegas, Gamma, rho);
}

/**
 * In the proposed example, the latent state x determines
 * the observational distribution:
 *  0: normal(mu, sigma)
 *  1: normal(-mu, sigma)
 */
double state_lpdf(double y, double abs_mu, double sigma, int state) {
  int x = state == 0 ? 1 : -1;
  double chi = (y - x * abs_mu) / sigma;
  return -0.5 * chi * chi - 0.5 * std::log(2 * M_PI) - std::log(sigma);
}

class hmm_test : public ::testing::Test {
 protected:
  void SetUp() override {
    n_states_ = 2;
    p1_init_ = 0.65;
    gamma1_ = 0.7;
    gamma2_ = 0.45;
    n_transitions_ = 10;
    abs_mu_ = 1;
    sigma_ = 1;

    Eigen::VectorXd rho(n_states_);
    rho << p1_init_, 1 - p1_init_;
    rho_ = rho;

    Eigen::MatrixXd Gamma(n_states_, n_states_);
    Gamma << gamma1_, 1 - gamma1_, gamma2_, 1 - gamma2_;
    Gamma_ = Gamma;

    Eigen::VectorXd obs_data(n_transitions_ + 1);
    obs_data << -0.3315914, -0.1655340, -0.7984021, 0.2364608, -0.4489722,
        2.1831438, -1.4778675, 0.8717423, -1.0370874, 0.1370296, 1.9786208;
    obs_data_ = obs_data;

    Eigen::MatrixXd log_omegas(n_states_, n_transitions_ + 1);
    for (int n = 0; n < n_transitions_ + 1; n++) {
      log_omegas.col(n)[0] = state_lpdf(obs_data[n], abs_mu_, sigma_, 0);
      log_omegas.col(n)[1] = state_lpdf(obs_data[n], abs_mu_, sigma_, 1);
    }
    log_omegas_ = log_omegas;
    log_omegas_zero_ = log_omegas.block(0, 0, n_states_, 1);

    std::vector<double> rho_unconstrained(n_states_ - 1);
    for (int i = 0; i < rho.size() - 1; i++)
      rho_unconstrained[i] = rho(i);
    rho_unconstrained_ = rho_unconstrained;

    Gamma_unconstrained_ = Gamma.block(0, 0, n_states_, n_states_ - 1);
  }

  int n_states_, n_transitions_;
  double abs_mu_, sigma_, p1_init_, gamma1_, gamma2_;

  Eigen::VectorXd rho_;
  Eigen::MatrixXd Gamma_;
  Eigen::VectorXd obs_data_;
  Eigen::MatrixXd log_omegas_;
  Eigen::MatrixXd log_omegas_zero_;

  // Construct "unconstrained" versions of rho and Gamma, without
  // the final element which can be determnied using the fact
  // the columns sum to 1. This allows us to do finite diff tests,
  // without violating the simplex constraint of rho and Gamma.
  std::vector<double> rho_unconstrained_;
  Eigen::MatrixXd Gamma_unconstrained_;
  stan::test::ad_tolerances tols_;
};
#endif
