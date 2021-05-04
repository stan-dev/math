#ifndef STAN_MATH_OPENCL_REV_OFFSET_MULTIPLIER_CONSTRAIN_HPP
#define STAN_MATH_OPENCL_REV_OFFSET_MULTIPLIER_CONSTRAIN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return the linearly transformed value for the specified unconstrained input
 * and specified offset and multiplier.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = mu + sigma * x\f$
 *
 * <p>where mu is the offset and sigma is the multiplier.
 *
 * <p>If the offset is zero and the multiplier is one this
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T type of unconstrained input
 * @tparam M type of offset
 * @tparam S type of multiplier
 * @param[in] A Unconstrained input
 * @param[in] mu offset of constrained output
 * @param[in] sigma multiplier of constrained output
 * @return linear transformed value corresponding to inputs
 * @throw std::domain_error if sigma <= 0
 * @throw std::domain_error if mu is not finite
 */
template <typename T, typename M, typename S,
          require_all_prim_or_rev_kernel_expression_t<T, M, S>* = nullptr,
          require_any_not_stan_scalar_t<T, M, S>* = nullptr,
          require_any_var_t<T, M, S>* = nullptr>
inline var_value<matrix_cl<double>> offset_multiplier_constrain(T&& A, M&& mu,
                                                                S&& sigma) {
  if (A.size() == 0) {
    return A;
  }
  arena_t<T> A_arena = std::forward<T>(A);
  arena_t<M> mu_arena = std::forward<M>(mu);
  arena_t<S> sigma_arena = std::forward<S>(sigma);
  return make_callback_var(
      offset_multiplier_constrain(value_of(A_arena), value_of(mu_arena),
                                  value_of(sigma_arena)),
      [A_arena, mu_arena,
       sigma_arena](vari_value<matrix_cl<double>>& res) mutable {
        adjoint_results(A_arena, mu_arena, sigma_arena) += expressions(
            elt_multiply(res.adj(), value_of(sigma_arena)), res.adj(),
            elt_multiply(res.adj(), value_of(A_arena)));
      });
}

/**
 * Return the linearly transformed value for the specified unconstrained input
 * and specified offset and multiplier.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = mu + sigma * x\f$
 *
 * <p>where mu is the offset and sigma is the multiplier.
 *
 * <p>If the offset is zero and the multiplier is one this
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T type of unconstrained input
 * @tparam M type of offset
 * @tparam S type of multiplier
 * @param[in] A Unconstrained input
 * @param[in] mu offset of constrained output
 * @param[in] sigma multiplier of constrained output
 * @param[in,out] lp Reference to log probability to increment.
 * @return linear transformed value corresponding to inputs
 * @throw std::domain_error if sigma <= 0
 * @throw std::domain_error if mu is not finite
 */
template <typename T, typename M, typename S,
          require_all_prim_or_rev_kernel_expression_t<T, M, S>* = nullptr,
          require_any_not_stan_scalar_t<T, M, S>* = nullptr,
          require_any_var_t<T, M, S>* = nullptr>
inline var_value<matrix_cl<double>> offset_multiplier_constrain(T&& A, M&& mu,
                                                                S&& sigma,
                                                                var& lp) {
  if (A.size() == 0) {
    return A;
  }
  arena_t<T> A_arena = std::forward<T>(A);
  arena_t<M> mu_arena = std::forward<M>(mu);
  arena_t<S> sigma_arena = std::forward<S>(sigma);

  double lp_inc = 0;
  auto res = offset_multiplier_constrain(value_of(A_arena), value_of(mu_arena),
                                         value_of(sigma_arena), lp_inc);
  lp += lp_inc;
  return make_callback_var(
      std::move(res), [A_arena, mu_arena, sigma_arena,
                       lp](vari_value<matrix_cl<double>>& res) mutable {
        adjoint_results(A_arena, mu_arena, sigma_arena) += expressions(
            elt_multiply(res.adj(), value_of(sigma_arena)), res.adj(),
            elt_multiply(res.adj(), value_of(A_arena))
                + elt_divide(lp.adj(), value_of(sigma_arena)));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
