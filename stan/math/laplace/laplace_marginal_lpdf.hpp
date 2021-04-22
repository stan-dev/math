#ifndef STAN_MATH_LAPLACE_LAPLACE_MARGINAL_LPDF_HPP
#define STAN_MATH_LAPLACE_LAPLACE_MARGINAL_LPDF_HPP

#include <stan/math/laplace/laplace_marginal.hpp>
#include <stan/math/laplace/laplace_likelihood_general.hpp>

namespace stan {
namespace math {
  /**
   * Wrapper function around the laplace_marginal function.
   * Returns the marginal density p(y | phi) by marginalizing out
   * the latent gaussian variable, with a Laplace approximation.
   * See the laplace_marginal function for more details.
   * The data y is assumed to be real.
   * The function is "overloaded" below for the int y and lpmf case.
   *
   * @tparam T0 The type of the initial guess, theta_0.
   * @tparam T1 The type for the global parameter, phi.
   * @tparam T2 The type of the auxiliary parameter, eta.
   * @tparam K The function which returns the prior covariance matrix.
   * @tparam F The function which returns the log likelihood.
   * @param[in] y fixed real data to be passed to the log likelihood.
   * @param[in] L_f a function which returns the log likelihood.
   * @param[in] eta non-marginalized parameters for the log likelihood.
   * @param[in] delta_int_f integer data to be passed to the log likelihood.
   * @param[in] K_f a function which returns the prior
   *              covariance for the marginalized out latent Gaussian.
   * @param[in] phi model parameters for the covariance function.
   * @param[in] x data for the covariance function.
   * @param[in] delta additional real data for the covariance matrix.
   * @param[in] delta_int_k additional int data for the covariance matrix.
   * @param[in] theta_0 initial guess for the Newton solver which returns
   *  the Laplace approximation.
   * @param[in] msgs_f message stream for the log likelihood function.
   * @param[in] msgs_k message stream for the covariance function.
   * @param[in] tolerance controls the convergence criterion when finding
   *            the mode in the Laplace approximation.
   * @param[in] max_num_steps maximum number of steps before the Newton solver
   *            breaks and returns an error.
   * @param[in] hessian_block_size the size of the block for a block-diagonal
   *              Hessian of the log likelihood. If 0, the Hessian is stored
   *              inside a vector. If the Hessian is dense, this should be the
   *              size of the Hessian.
   * @param[in] compute_W_root if 1, the Newton solver computes the root of W,
   *              the negative Hessian of the log likelihood, which leads to
   *              efficient computation. Else, a more general but slower solver
   *              is used.
   */
  /*
  template <bool propto, typename T0, typename T1, typename T2,
            typename Tx, typename K, typename L>
  stan::return_type_t<T1, T2> laplace_marginal_lpdf
    (const Eigen::VectorXd& y,
     const L& L_f,
     const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta,
     const std::vector<int>& delta_int_L,
     const K& K_f,
     const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
     const Tx& x,
     const std::vector<double>& delta_K,
     const std::vector<int>& delta_int_K,
     const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
     std::ostream* msgs_L = nullptr,
     std::ostream* msgs_K = nullptr,
     double tolerance = 1e-6,
     long int max_num_steps = 100,
     int hessian_block_size = 0,
     int compute_W_root = 1) {

    return laplace_marginal_density(
      diff_likelihood<L>(L_f, y, delta_int_L, msgs_L),
      K_f, phi, eta, x, delta_K, delta_int_K,
      theta_0, msgs_K, tolerance, max_num_steps,
      hessian_block_size, compute_W_root);
  }  */

  template <bool propto, typename T0, typename T1, typename T2,
            typename Tx, typename K, typename L>
  stan::return_type_t<T1, T2> laplace_marginal_lpdf
    (const Eigen::VectorXd& y,
     const L& L_f,
     const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta,
     const std::vector<int>& delta_int_L,
     const K& K_f,
     const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
     const Tx& x,
     const std::vector<double>& delta_K,
     const std::vector<int>& delta_int_K,
     const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
     double tolerance = 1e-6,
     long int max_num_steps = 100,
     int hessian_block_size = 0,
     int solver = 1,
     int do_line_search = 1,
     int max_steps_line_search = 10,
     std::ostream* msgs = nullptr) {
    // TEST: provisional signature to agree with parser.

    return laplace_marginal_density(
      diff_likelihood<L>(L_f, y, delta_int_L, msgs),
      K_f, phi, eta, x, delta_K, delta_int_K,
      theta_0, msgs, tolerance, max_num_steps,
      hessian_block_size, solver,
      do_line_search, max_steps_line_search);
  }

  /**
   * Overloaded function for lpmf case. The first argument
   * is now a std::vector of interger and an Eigen::VectorXd
   * of double is passed as data.
   */
   template <bool propto, typename T0, typename T1, typename T2,
             typename Tx, typename K, typename L>
   stan::return_type_t<T1, T2> laplace_marginal_lpmf
     (const std::vector<int>& y,
      const L& L_f,
      const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta,
      const Eigen::VectorXd& delta_L,
      const K& K_f,
      const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
      const Tx& x,
      const std::vector<double>& delta_K,
      const std::vector<int>& delta_int_K,
      const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
      double tolerance = 1e-6,
      long int max_num_steps = 100,
      int hessian_block_size = 0,
      int solver = 1,
      int do_line_search = 1,
      int max_steps_line_search = 10,
      std::ostream* msgs = nullptr) {

    return laplace_marginal_lpdf<propto>(delta_L, L_f, eta, y,
                                 K_f, phi, x, delta_K, delta_int_K,
                                 theta_0, tolerance,
                                 max_num_steps,
                                 hessian_block_size,
                                 solver, do_line_search,
                                 max_steps_line_search, msgs);
  }
}  // namespace math
}  // namespace stan

#endif
