#ifndef STAN_MATH_PRIM_PROB_LOGLOGISTIC_CDF_HPP
#define STAN_MATH_PRIM_PROB_LOGLOGISTIC_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/F32.hpp>
#include <stan/math/prim/fun/grad_F32.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Docs describing the templates and arguments
 */
template <typename T_y, typename T_alpha, typename T_beta>
return_type_t<T_y, T_alpha, T_beta> loglogistic_cdf(const T_y& y,
    const T_alpha& alpha,
    const T_beta& beta) {
  using T_partials_return = partials_return_t<T_y, T_alpha, T_beta>;
  using std::exp;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_alpha>;
  using T_beta_ref = ref_type_t<T_beta>;
  constexpr const char* function = "loglogistic_cdf";
  check_consistent_sizes(function, "y", y, "Shape", alpha, "Scale", beta);
  if (size_zero(n, N, alpha, beta)) {
    return 1.0;
  }

  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_nonnegative(function, "y", y);
  check_nonnegative(function, "Shape", alpha_ref);
  check_nonnegative(function, "Scale", beta_ref);

  T_partials_return logp(1.0);
  operands_and_partials<T_y_ref, T_alpha_ref, T_beta_ref> ops_partials(y_ref, alpha_ref, beta_ref);

  scalar_seq_view<T_y_ref> N_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t max_size_seq_view = max_size(y, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(y); i++) {
    if (y_vec.val(i) == NEGATIVE_INFTY) {
      return ops_partials.build(0.0);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; ++i) {
    auto y_val = y_vec.val(i);
    auto alpha_val = alpha_vec.val(i);
    auto beta_val = beta_vec.val(i);
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[i] += ;// Partial for y;
    }
    if (!is_constant_all<T_alpha>::value) {
      ops_partials.edge1_.partials_[i] += ;// Partial for alpha;
    }
    if (!is_constant_all<T_beta>::value) {
      ops_partials.edge2_.partials_[i] += ;// Partial for beta;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
