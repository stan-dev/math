#ifndef STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_LPDF_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_LPDF_HPP

#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/laplace_marginal_density.hpp>

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
 * @tparam LFun The function which returns the log likelihood.
 * @tparam LArgs A tuple of arguments to the log likelihood.
 * @tparam EtaVec The type of the auxiliary parameter, eta.
 * @tparam CovarFun The function which returns the prior covariance matrix.
 * @tparam Theta0 The type of the initial guess, theta_0.
 * @tparam Args Arguments supplied to the CovarFun
 * @param[in] L_f a function which returns the log likelihood.
 * @param[in] l_args A tuple of arguments to pass to the log likelihood
 * @param[in] eta non-marginalized parameters for the log likelihood.
 * @param[in] theta_0 initial guess for the Newton solver which returns
 *  the Laplace approximation.
 * @param[in] K_f a function which returns the prior
 *              covariance for the marginalized out latent Gaussian.
 * @param[in] tolerance controls the convergence criterion when finding
 *            the mode in the Laplace approximation.
 * @param[in] max_num_steps maximum number of steps before the Newton solver
 *            breaks and returns an error.
 * @param[in] hessian_block_size the size of the block for a block-diagonal
 *              Hessian of the log likelihood. If 0, the Hessian is stored
 *              inside a vector. If the Hessian is dense, this should be the
 *              size of the Hessian.
 * @param[in] solver
 * @param[in] max_steps_line_search
 * @param[in] msgs message stream for the covariance and likelihood function.
 * @param[in] args A parameter pack of arguments to be send to the covariance
 * function
 */
template <bool propto = false, typename LFun, typename LArgs, typename EtaVec,
          typename CovarFun, typename Theta0, typename... Args,
          require_all_eigen_vector_t<EtaVec, Theta0>* = nullptr>
inline auto laplace_marginal_tol_lpdf(
    LFun&& L_f, LArgs&& l_args, const EtaVec& eta, const Theta0& theta_0,
    CovarFun&& K_f, double tolerance, int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs, Args&&... args) {
  // TEST: provisional signature to agree with parser.
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_marginal_density(std::forward<LFun>(L_f),
                                  std::forward<LArgs>(l_args),
                                  std::forward<CovarFun>(K_f), eta, theta_0,
                                  msgs, ops, std::forward<Args>(args)...);
}

// no eta
template <bool propto = false, typename LArgs, typename LFun, typename CovarFun,
          typename Theta0, typename... Args,
          require_all_eigen_vector_t<Theta0>* = nullptr>
inline auto laplace_marginal_tol_lpdf(LFun&& L_f, LArgs&& l_args,
                                      const Theta0& theta_0, CovarFun&& K_f,
                                      double tolerance, int64_t max_num_steps,
                                      const int hessian_block_size,
                                      const int solver,
                                      const int max_steps_line_search,
                                      std::ostream* msgs, Args&&... args) {
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  Eigen::Matrix<double, 0, 0> eta;
  return laplace_marginal_density(std::forward<LFun>(L_f),
                                  std::forward<LArgs>(l_args),
                                  std::forward<CovarFun>(K_f), eta, theta_0,
                                  msgs, ops, std::forward<Args>(args)...);
}

/**
 * Overloaded function for lpmf case. The first argument
 * is now a std::vector of interger and an Eigen::VectorXd
 * of double is passed as data.
 */
template <bool propto = false, typename LFun, typename LArgs, typename EtaVec,
          typename CovarFun, typename Theta0, typename... Args>
inline auto laplace_marginal_tol_lpmf(
    LFun&& L_f, LArgs&& l_args, const EtaVec& eta, const Theta0& theta_0,
    CovarFun&& K_f, const double tolerance, const int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs, Args&&... args) {
  return laplace_marginal_tol_lpdf<propto>(
      std::forward<LFun>(L_f), std::forward<LArgs>(l_args), eta, tolerance,
      max_num_steps, hessian_block_size, solver, max_steps_line_search, theta_0,
      std::forward<CovarFun>(K_f), msgs, std::forward<Args>(args)...);
}

template <bool propto = false, typename LFun, typename LArgs, typename EtaVec,
          typename CovarFun, typename Theta0, typename... Args>
inline auto laplace_marginal_tol_lpmf(
    LFun&& L_f, LArgs&& l_args, const Theta0& theta_0, CovarFun&& K_f,
    const double tolerance, const int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs, Args&&... args) {
  return laplace_marginal_tol_lpdf<propto>(
      std::forward<LFun>(L_f), std::forward<LArgs>(l_args), theta_0,
      std::forward<CovarFun>(K_f), tolerance, max_num_steps, hessian_block_size,
      solver, max_steps_line_search, msgs, std::forward<Args>(args)...);
}

/**
 * Overloaded function for lpmf case. The first argument
 * is now a std::vector of interger and an Eigen::VectorXd
 * of double is passed as data.
 */
template <bool propto = false, typename LFun, typename LArgs, typename CovarFun,
          typename Theta0, typename... Args>
inline auto laplace_marginal_tol_lpmf(
    LFun&& L_f, LArgs&& l_args, const Theta0& theta_0, CovarFun&& K_f,
    const double tolerance, const int64_t max_num_steps,
    const int hessian_block_size, const int solver,
    const int max_steps_line_search, std::ostream* msgs, Args&&... args) {
  return laplace_marginal_tol_lpdf<propto>(
      std::forward<LFun>(L_f), std::forward<LArgs>(l_args), theta_0,
      std::forward<CovarFun>(K_f), tolerance, max_num_steps, hessian_block_size,
      solver, max_steps_line_search, msgs, std::forward<Args>(args)...);
}

template <bool propto = false, typename LFun, typename LArgs, typename EtaVec,
          typename CovarFun, typename Theta0, typename... Args,
          require_all_eigen_vector_t<EtaVec, Theta0>* = nullptr>
inline auto laplace_marginal_lpdf(LFun&& L_f, LArgs&& l_args, const EtaVec& eta,
                                  const Theta0& theta_0, CovarFun&& K_f,
                                  std::ostream* msgs, Args&&... args) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(std::forward<LFun>(L_f),
                                  std::forward<LArgs>(l_args),
                                  std::forward<CovarFun>(K_f), eta, theta_0,
                                  msgs, ops, std::forward<Args>(args)...);
}

template <bool propto = false, typename LFun, typename LArgs, typename CovarFun,
          typename Theta0, typename... Args,
          require_all_eigen_vector_t<Theta0>* = nullptr>
inline auto laplace_marginal_lpdf(LFun&& L_f, LArgs&& l_args,
                                  const Theta0& theta_0, CovarFun&& K_f,
                                  std::ostream* msgs, Args&&... args) {
  // TEST: provisional signature to agree with parser.
  Eigen::Matrix<double, 0, 0> eta;
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(std::forward<LFun>(L_f),
                                  std::forward<LArgs>(l_args),
                                  std::forward<CovarFun>(K_f), eta, theta_0,
                                  msgs, ops, std::forward<Args>(args)...);
}

template <bool propto = false, typename LFun, typename LArgs, typename EtaVec,
          typename CovarFun, typename Theta0, typename... Args,
          require_all_eigen_vector_t<EtaVec, Theta0>* = nullptr>
inline auto laplace_marginal_lpmf(LFun&& L_f, LArgs&& l_args, const EtaVec& eta,
                                  const Theta0& theta_0, CovarFun&& K_f,
                                  std::ostream* msgs, Args&&... args) {
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(std::forward<LFun>(L_f),
                                  std::forward<LArgs>(l_args),
                                  std::forward<CovarFun>(K_f), eta, theta_0,
                                  msgs, ops, std::forward<Args>(args)...);
}

template <bool propto = false, typename LFun, typename LArgs, typename CovarFun,
          typename Theta0, typename... Args,
          require_all_eigen_vector_t<Theta0>* = nullptr>
inline auto laplace_marginal_lpmf(LFun&& L_f, LArgs&& l_args,
                                  const Theta0& theta_0, CovarFun&& K_f,
                                  std::ostream* msgs, Args&&... args) {
  // TEST: provisional signature to agree with parser.
  Eigen::Matrix<double, 0, 0> eta;
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(std::forward<LFun>(L_f),
                                  std::forward<LArgs>(l_args),
                                  std::forward<CovarFun>(K_f), eta, theta_0,
                                  msgs, ops, std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
