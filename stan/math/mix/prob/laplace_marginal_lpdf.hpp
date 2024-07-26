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
template <bool propto = false, typename YVec, typename LFun, typename EtaVec,
          typename CovarFun, typename Theta0, typename... Args,
          require_all_eigen_vector_t<YVec, EtaVec, Theta0>* = nullptr>
inline auto laplace_marginal_tol_lpdf(
    const YVec& y, LFun&& L_f, const EtaVec& eta,
    const std::vector<int>& delta_int_L, double tolerance,
    long int max_num_steps, const int hessian_block_size, const int solver,
    const int max_steps_line_search, const Theta0& theta_0, CovarFun&& K_f,
    std::ostream* msgs, Args&&... args) {
  // TEST: provisional signature to agree with parser.
    laplace_options ops{hessian_block_size, solver,
    max_steps_line_search, tolerance, max_num_steps};
  return laplace_marginal_density(
      laplace_likelihood<LFun>(std::forward<LFun>(L_f), y, delta_int_L, msgs),
      std::forward<CovarFun>(K_f), eta, theta_0, msgs, ops,
      std::forward<Args>(args)...);
}

template <bool propto = false, typename YVec, typename LFun,
          typename CovarFun, typename Theta0, typename... Args,
          require_all_eigen_vector_t<YVec, Theta0>* = nullptr>
inline auto laplace_marginal_tol_lpdf(
    const YVec& y, LFun&& L_f,
    const std::vector<int>& delta_int_L, double tolerance,
    long int max_num_steps, const int hessian_block_size, const int solver,
    const int max_steps_line_search, const Theta0& theta_0, CovarFun&& K_f,
    std::ostream* msgs, Args&&... args) {
    laplace_options ops{hessian_block_size, solver,
    max_steps_line_search, tolerance, max_num_steps};
    Eigen::Matrix<double, 0, 0> eta;
  return laplace_marginal_density(
      laplace_likelihood<LFun>(std::forward<LFun>(L_f), y, delta_int_L, msgs),
      std::forward<CovarFun>(K_f), eta, theta_0, msgs, ops,
      std::forward<Args>(args)...);
}

/**
 * Overloaded function for lpmf case. The first argument
 * is now a std::vector of interger and an Eigen::VectorXd
 * of double is passed as data.
 */
template <bool propto = false, typename LFun, typename EtaVec,
          typename CovarFun, typename DeltaLVec, typename Theta0,
          typename... Args>
inline auto laplace_marginal_tol_lpmf(
    const std::vector<int>& y, LFun&& L_f, const EtaVec& eta,
    const DeltaLVec& delta_L, const double tolerance,
    const long int max_num_steps, const int hessian_block_size,
    const int solver, const int max_steps_line_search, const Theta0& theta_0,
    CovarFun&& K_f, std::ostream* msgs, Args&&... args) {
  return laplace_marginal_tol_lpdf<propto>(
      delta_L, std::forward<LFun>(L_f), eta, y, tolerance, max_num_steps,
      hessian_block_size, solver, max_steps_line_search, theta_0,
      std::forward<CovarFun>(K_f), msgs, std::forward<Args>(args)...);
}

/**
 * Overloaded function for lpmf case. The first argument
 * is now a std::vector of interger and an Eigen::VectorXd
 * of double is passed as data.
 */
template <bool propto = false, typename LFun,
          typename CovarFun, typename DeltaLVec, typename Theta0,
          typename... Args>
inline auto laplace_marginal_tol_lpmf(
    const std::vector<int>& y, LFun&& L_f,
    const DeltaLVec& delta_L, const double tolerance,
    const long int max_num_steps, const int hessian_block_size,
    const int solver, const int max_steps_line_search, const Theta0& theta_0,
    CovarFun&& K_f, std::ostream* msgs, Args&&... args) {
    Eigen::Matrix<double, 0, 0> eta;
  return laplace_marginal_tol_lpdf<propto>(
      delta_L, std::forward<LFun>(L_f), eta, y, tolerance, max_num_steps,
      hessian_block_size, solver, max_steps_line_search, theta_0,
      std::forward<CovarFun>(K_f), msgs, std::forward<Args>(args)...);
}

template <bool propto = false, typename YVec, typename LFun, typename EtaVec,
          typename CovarFun, typename Theta0, typename... Args,
          require_all_eigen_vector_t<YVec, EtaVec, Theta0>* = nullptr>
inline auto laplace_marginal_lpdf(const YVec& y, LFun&& L_f, const EtaVec& eta,
                                  const std::vector<int>& delta_int_L,
                                  const Theta0& theta_0, CovarFun&& K_f,
                                  std::ostream* msgs, Args&&... args) {
  // TEST: provisional signature to agree with parser.
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(
      laplace_likelihood<LFun>(std::forward<LFun>(L_f), y, delta_int_L, msgs),
      std::forward<CovarFun>(K_f), eta, theta_0, msgs, ops,
      std::forward<Args>(args)...);
}

template <bool propto = false, typename YVec, typename LFun,
          typename CovarFun, typename Theta0, typename... Args,
          require_all_eigen_vector_t<YVec, Theta0>* = nullptr>
inline auto laplace_marginal_lpdf(const YVec& y, LFun&& L_f,
                                  const std::vector<int>& delta_int_L,
                                  const Theta0& theta_0, CovarFun&& K_f,
                                  std::ostream* msgs, Args&&... args) {
  // TEST: provisional signature to agree with parser.
    Eigen::Matrix<double, 0, 0> eta;
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(
      laplace_likelihood<LFun>(std::forward<LFun>(L_f), y, delta_int_L, msgs),
      std::forward<CovarFun>(K_f), eta, theta_0, msgs, ops,
      std::forward<Args>(args)...);
}

/**
 * Overloaded function for lpmf case. The first argument
 * is now a std::vector of interger and an Eigen::VectorXd
 * of double is passed as data.
 */
template <bool propto = false, typename LFun, typename EtaVec,
          typename CovarFun, typename DeltaLVec, typename Theta0,
          typename... Args>
inline auto laplace_marginal_lpmf(const std::vector<int>& y, LFun&& L_f,
                                  const EtaVec& eta, const DeltaLVec& delta_L,
                                  const Theta0& theta_0, CovarFun&& K_f,
                                  std::ostream* msgs, Args&&... args) {
  constexpr double tolerance = 1e-6;
  constexpr long int max_num_steps = 100;
  constexpr int hessian_block_size = 1;
  constexpr int solver = 1;
  constexpr int max_steps_line_search = 0;
  return laplace_marginal_tol_lpdf<propto>(
      delta_L, std::forward<LFun>(L_f), eta, y, tolerance, max_num_steps,
      hessian_block_size, solver, max_steps_line_search, theta_0,
      std::forward<CovarFun>(K_f), msgs, std::forward<Args>(args)...);
}

/**
 * Overloaded function for lpmf case. The first argument
 * is now a std::vector of interger and an Eigen::VectorXd
 * of double is passed as data.
 */
template <bool propto = false, typename LFun,
          typename CovarFun, typename DeltaLVec, typename Theta0,
          typename... Args>
inline auto laplace_marginal_lpmf(const std::vector<int>& y, LFun&& L_f,
                                  const DeltaLVec& delta_L,
                                  const Theta0& theta_0, CovarFun&& K_f,
                                  std::ostream* msgs, Args&&... args) {
  constexpr double tolerance = 1e-6;
  constexpr long int max_num_steps = 100;
  constexpr int hessian_block_size = 1;
  constexpr int solver = 1;
  constexpr int max_steps_line_search = 0;
  Eigen::Matrix<double, 0, 0> eta;
  return laplace_marginal_tol_lpdf<propto>(
      delta_L, std::forward<LFun>(L_f), eta, y, tolerance, max_num_steps,
      hessian_block_size, solver, max_steps_line_search, theta_0,
      std::forward<CovarFun>(K_f), msgs, std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
