#ifndef STAN_MATH_PRIM_MAT_PROB_WIENER_LOG_HPP
#define STAN_MATH_PRIM_MAT_PROB_WIENER_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/wiener_lpdf.hpp>
#include <boost/math/tools/promotion.hpp>

namespace stan {
namespace math {

/**
 * The log of the first passage time density function for a (Wiener)
 *  drift diffusion model for the given \f$y\f$,
 * boundary separation \f$\alpha\f$, nondecision time \f$\tau\f$,
 * relative bias \f$\beta\f$, and drift rate \f$\delta\f$.
 * \f$\alpha\f$ and \f$\tau\f$ must be greater than 0, and
 * \f$\beta\f$ must be between 0 and 1. \f$y\f$ should contain
 * reaction times in seconds (strictly positive) with
 * upper-boundary responses.
 *
 * @deprecated use <code>wiener_lpdf</code>
 *
 * @param y A scalar variate.
 * @param alpha The boundary separation.
 * @param tau The nondecision time.
 * @param beta The relative bias.
 * @param delta The drift rate.
 * @return The log of the Wiener first passage time density of
 *  the specified arguments.
 */
template <bool propto, typename T_y, typename T_alpha, typename T_tau,
          typename T_beta, typename T_delta>
inline auto wiener_log(T_y&& y, T_alpha&& alpha, T_tau&& tau,
                       T_beta&& beta, T_delta&& delta) {
  return wiener_lpdf<propto>(std::forward<T_y>(y), std::forward<T_alpha>(alpha),
   std::forward<T_tau>(tau), std::forward<T_beta>(beta), std::forward<T_delta>(delta));
}

/**
 * @deprecated use <code>wiener_lpdf</code>
 */
template <typename T_y, typename T_alpha, typename T_tau, typename T_beta,
          typename T_delta>
inline auto wiener_log(T_y&& y, T_alpha&& alpha, T_tau&& tau,
                       T_beta&& beta, T_delta&& delta) {
  return wiener_lpdf(std::forward<T_y>(y), std::forward<T_alpha>(alpha),
   std::forward<T_tau>(tau), std::forward<T_beta>(beta), std::forward<T_delta>(delta));
}

}  // namespace math
}  // namespace stan
#endif
