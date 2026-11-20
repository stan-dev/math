#ifndef STAN_MATH_PRIM_PROB_CHI_SQUARE_CDF_HPP
#define STAN_MATH_PRIM_PROB_CHI_SQUARE_CDF_HPP

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
 * The cumulative distribution function of the chi-square distribution.
 * It is computed from the log of the regularized lower
 * incomplete gamma function, so it stays finite where the probability
 * underflows.
 *
 * @tparam T_y type of y
 * @tparam T_dof type of nu
 * @param y random variable
 * @param nu degrees of freedom parameter
 * @return probability, or product of probabilities for containers
 * @throw std::domain_error if a parameter is out of its support
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_y, typename T_dof>
inline return_type_t<T_y, T_dof> chi_square_cdf(const T_y& y, const T_dof& nu) {
  using T_partials_return = partials_return_t<T_y, T_dof>;
  using T_y_ref = ref_type_t<T_y>;
  using T_dof_ref = ref_type_t<T_dof>;
  static constexpr const char* function = "chi_square_cdf";
  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu);
  T_y_ref y_ref = y;
  T_dof_ref nu_ref = nu;
  check_not_nan(function, "Random variable", y_ref);
  check_nonnegative(function, "Random variable", y_ref);
  check_positive_finite(function, "Degrees of freedom parameter", nu_ref);

  if (size_zero(y, nu)) {
    return 1.0;
  }

  T_partials_return log_P(0.0);
  auto ops_partials = make_partials_propagator(y_ref, nu_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_dof_ref> nu_vec(nu_ref);
  const size_t N = max_size(y, nu);

  // An element with probability 0 (z = 0) makes the
  // result 0; the gradients are treated as zero
  for (size_t n = 0; n < N; n++) {
    const double y_rec = value_of_rec(y_vec.val(n));
    if (0.5 * y_rec == 0) {
      return ops_partials.build(0.0);
    }
  }

  constexpr bool need_d_a = is_autodiff_v<T_dof>;
  constexpr bool need_d_z = is_autodiff_v<T_y>;
  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_val = y_vec.val(n);
    const T_partials_return nu_val = nu_vec.val(n);
    const T_partials_return z = 0.5 * y_val;
    const double z_rec = value_of_rec(z);
    // An element with probability 1 contributes nothing
    if (z_rec == INFTY) {
      continue;
    }
    const auto tail
        = internal::log_gamma_tail<false, need_d_a, need_d_z>(0.5 * nu_val, z);
    log_P += tail.log_prob;
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials)[n] += 0.5 * tail.d_z;
    }
    if constexpr (is_autodiff_v<T_dof>) {
      partials<1>(ops_partials)[n] += 0.5 * tail.d_a;
    }
  }
  const T_partials_return P = exp(log_P);
  if constexpr (is_autodiff_v<T_y>) {
    for (size_t n = 0; n < stan::math::size(y); ++n) {
      partials<0>(ops_partials)[n] *= P;
    }
  }
  if constexpr (is_autodiff_v<T_dof>) {
    for (size_t n = 0; n < stan::math::size(nu); ++n) {
      partials<1>(ops_partials)[n] *= P;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
