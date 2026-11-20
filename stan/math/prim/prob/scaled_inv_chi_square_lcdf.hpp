#ifndef STAN_MATH_PRIM_PROB_SCALED_INV_CHI_SQUARE_LCDF_HPP
#define STAN_MATH_PRIM_PROB_SCALED_INV_CHI_SQUARE_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log_gamma_q_dgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log cumulative distribution function of the scaled inverse
 * chi-square distribution.
 * It is computed from the log of the regularized upper
 * incomplete gamma function, so it stays finite where the probability
 * underflows.
 *
 * @tparam T_y type of y
 * @tparam T_dof type of nu
 * @tparam T_scale type of s
 * @param y random variable
 * @param nu degrees of freedom parameter
 * @param s scale parameter
 * @return log probability, or sum of log probabilities for containers
 * @throw std::domain_error if a parameter is out of its support
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_y, typename T_dof, typename T_scale>
inline return_type_t<T_y, T_dof, T_scale> scaled_inv_chi_square_lcdf(
    const T_y& y, const T_dof& nu, const T_scale& s) {
  using T_partials_return = partials_return_t<T_y, T_dof, T_scale>;
  using T_y_ref = ref_type_t<T_y>;
  using T_dof_ref = ref_type_t<T_dof>;
  using T_scale_ref = ref_type_t<T_scale>;
  static constexpr const char* function = "scaled_inv_chi_square_lcdf";
  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu, "Scale parameter",
                         s);
  T_y_ref y_ref = y;
  T_dof_ref nu_ref = nu;
  T_scale_ref s_ref = s;
  check_nonnegative(function, "Random variable", y_ref);
  check_positive_finite(function, "Degrees of freedom parameter", nu_ref);
  check_positive_finite(function, "Scale parameter", s_ref);

  if (size_zero(y, nu, s)) {
    return 0;
  }

  T_partials_return log_P(0.0);
  auto ops_partials = make_partials_propagator(y_ref, nu_ref, s_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_dof_ref> nu_vec(nu_ref);
  scalar_seq_view<T_scale_ref> s_vec(s_ref);
  const size_t N = max_size(y, nu, s);

  // An element with probability 0 (z = inf) makes the
  // result -inf; the gradients are treated as zero
  for (size_t n = 0; n < N; n++) {
    const double y_rec = value_of_rec(y_vec.val(n));
    const double nu_rec = value_of_rec(nu_vec.val(n));
    const double s_rec = value_of_rec(s_vec.val(n));
    if (0.5 * nu_rec * s_rec * s_rec / y_rec == INFTY) {
      return ops_partials.build(negative_infinity());
    }
  }

  constexpr bool need_d_a = is_autodiff_v<T_dof>;
  constexpr bool need_d_z = is_any_autodiff_v<T_y, T_dof, T_scale>;
  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_val = y_vec.val(n);
    const T_partials_return nu_val = nu_vec.val(n);
    const T_partials_return s_val = s_vec.val(n);
    const T_partials_return z = 0.5 * nu_val * s_val * s_val / y_val;
    const double z_rec = value_of_rec(z);
    // An element with probability 1 contributes nothing
    if (z_rec == 0) {
      continue;
    }
    const auto tail
        = internal::log_gamma_tail<true, need_d_a, need_d_z>(0.5 * nu_val, z);
    log_P += tail.log_prob;
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials)[n] += -z / y_val * tail.d_z;
    }
    if constexpr (is_autodiff_v<T_dof>) {
      partials<1>(ops_partials)[n]
          += 0.5 * tail.d_a + 0.5 * s_val * s_val / y_val * tail.d_z;
    }
    if constexpr (is_autodiff_v<T_scale>) {
      partials<2>(ops_partials)[n] += nu_val * s_val / y_val * tail.d_z;
    }
  }
  return ops_partials.build(log_P);
}

}  // namespace math
}  // namespace stan
#endif
