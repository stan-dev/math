#ifndef STAN_MATH_PRIM_PROB_PARETO_CDF_HPP
#define STAN_MATH_PRIM_PROB_PARETO_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_scale, typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_scale, T_shape>* = nullptr>
return_type_t<T_y, T_scale, T_shape> pareto_cdf(const T_y& y,
                                                const T_scale& y_min,
                                                const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_scale, T_shape>;
  using T_y_ref = ref_type_t<T_y>;
  using T_y_min_ref = ref_type_t<T_scale>;
  using T_alpha_ref = ref_type_t<T_shape>;
  using std::exp;
  using std::log;
  static const char* function = "pareto_cdf";
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         y_min, "Shape parameter", alpha);
  if (size_zero(y, y_min, alpha)) {
    return 1.0;
  }

  T_y_ref y_ref = y;
  T_y_min_ref y_min_ref = y_min;
  T_alpha_ref alpha_ref = alpha;

  check_nonnegative(function, "Random variable", y_ref);
  check_positive_finite(function, "Scale parameter", y_min_ref);
  check_positive_finite(function, "Shape parameter", alpha_ref);

  T_partials_return P(1.0);
  auto ops_partials = make_partials_propagator(y_ref, y_min_ref, alpha_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_y_min_ref> y_min_vec(y_min_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  size_t N = max_size(y, y_min, alpha);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(y); i++) {
    if (y_vec.val(i) < y_min_vec.val(i)) {
      return ops_partials.build(0.0);
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (y_vec.val(n) == INFTY) {
      continue;
    }

    const T_partials_return log_dbl = log(y_min_vec.val(n) / y_vec.val(n));
    const T_partials_return y_min_inv_dbl = 1.0 / y_min_vec.val(n);
    const T_partials_return alpha_dbl = alpha_vec.val(n);

    const T_partials_return Pn = 1.0 - exp(alpha_dbl * log_dbl);

    P *= Pn;

    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials)[n]
          += alpha_dbl * y_min_inv_dbl * exp((alpha_dbl + 1) * log_dbl) / Pn;
    }
    if (!is_constant_all<T_scale>::value) {
      partials<1>(ops_partials)[n]
          += -alpha_dbl * y_min_inv_dbl * exp(alpha_dbl * log_dbl) / Pn;
    }
    if (!is_constant_all<T_shape>::value) {
      partials<2>(ops_partials)[n] += -exp(alpha_dbl * log_dbl) * log_dbl / Pn;
    }
  }

  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < stan::math::size(y); ++n) {
      partials<0>(ops_partials)[n] *= P;
    }
  }
  if (!is_constant_all<T_scale>::value) {
    for (size_t n = 0; n < stan::math::size(y_min); ++n) {
      partials<1>(ops_partials)[n] *= P;
    }
  }
  if (!is_constant_all<T_shape>::value) {
    for (size_t n = 0; n < stan::math::size(alpha); ++n) {
      partials<2>(ops_partials)[n] *= P;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
