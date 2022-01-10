#ifndef STAN_MATH_REV_PROB_EXPONENTIAL_QF_HPP
#define STAN_MATH_REV_PROB_EXPONENTIAL_QF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/prob/exponential_qf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The quantile of an exponential density for p with the specified
 * inverse scale parameter.
 * Inverse scale parameter must be greater than 0.
 * p must be bounded by 0 and 1.
 *
 \f[
   \frac{\partial }{\partial p} = \frac{1}{\beta - \beta p} \\
   \frac{\partial }{\partial \beta} = \frac{\log(1 - p)}{\beta^2}
 \f]
 *
 * @tparam Tp type of probability input
 * @tparam Tbeta type of inverse scale
 * @param p A scalar variable.
 * @param beta Inverse scale parameter.
 * @throw std::domain_error if beta is not greater than 0.
 * @throw std::domain_error if y is not greater than or equal to 0.
 */
template <typename Tp, typename Tbeta, require_any_st_var<Tp, Tbeta>* = nullptr>
inline auto exponential_qf(const Tp& p, const Tbeta& beta) {
  using inner_ret_type = decltype(exponential_qf(value_of(p), value_of(beta)));
  using ret_type = return_var_matrix_t<inner_ret_type, Tp, Tbeta>;

  arena_t<Tp> arena_p = p;
  arena_t<Tbeta> arena_beta = beta;
  arena_t<ret_type> ret
      = exponential_qf(value_of(arena_p), value_of(arena_beta));
  reverse_pass_callback([ret, arena_p, arena_beta]() mutable {
    decltype(auto) p_val = as_array_or_scalar(value_of(arena_p));
    decltype(auto) beta_val = as_array_or_scalar(value_of(arena_beta));
    decltype(auto) ret_adj = as_array_or_scalar(ret.adj());
    if (!is_constant<Tp>::value) {
      as_array_or_scalar(forward_as<promote_scalar_t<var, Tp>>(arena_p).adj())
          += ret_adj * inv(beta_val - beta_val * p_val);
    }
    if (!is_constant<Tbeta>::value) {
      as_array_or_scalar(
          forward_as<promote_scalar_t<var, Tbeta>>(arena_beta).adj())
          += ret_adj * log1m(p_val) / square(beta_val);
    }
  });
  return ret_type(ret);
}

}  // namespace math
}  // namespace stan
#endif
