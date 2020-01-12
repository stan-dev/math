#ifndef STAN_MATH_PRIM_PROB_PARETO_LPDF_HPP
#define STAN_MATH_PRIM_PROB_PARETO_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

// Pareto(y|y_m, alpha)  [y > y_m;  y_m > 0;  alpha > 0]
template <bool propto, typename T_y, typename T_scale, typename T_shape>
return_type_t<T_y, T_scale, T_shape> pareto_lpdf(const T_y& y,
                                                 const T_scale& y_min,
                                                 const T_shape& alpha) {
  static const char* function = "pareto_lpdf";
  using T_partials_return = partials_return_t<T_y, T_scale, T_shape>;
  using std::log;
  check_not_nan(function, "Random variable", y);
  check_positive_finite(function, "Scale parameter", y_min);
  check_positive_finite(function, "Shape parameter", alpha);
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         y_min, "Shape parameter", alpha);
  if (size_zero(y, y_min, alpha)) {
    return 0;
  }

  if (!include_summand<propto, T_y, T_scale, T_shape>::value) {
    return 0;
  }

  T_partials_return logp(0);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_scale> y_min_vec(y_min);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  size_t N = max_size(y, y_min, alpha);

  for (size_t n = 0; n < N; n++) {
    if (y_vec[n] < y_min_vec[n]) {
      return LOG_ZERO;
    }
  }

  operands_and_partials<T_y, T_scale, T_shape> ops_partials(y, y_min, alpha);

  VectorBuilder<include_summand<propto, T_y, T_shape>::value, T_partials_return,
                T_y>
      log_y(size(y));
  if (include_summand<propto, T_y, T_shape>::value) {
    for (size_t n = 0; n < size(y); n++) {
      log_y[n] = log(value_of(y_vec[n]));
    }
  }

  VectorBuilder<!is_constant_all<T_y, T_shape>::value, T_partials_return, T_y>
      inv_y(size(y));
  if (!is_constant_all<T_y, T_shape>::value) {
    for (size_t n = 0; n < size(y); n++) {
      inv_y[n] = 1 / value_of(y_vec[n]);
    }
  }

  VectorBuilder<include_summand<propto, T_scale, T_shape>::value,
                T_partials_return, T_scale>
      log_y_min(size(y_min));
  if (include_summand<propto, T_scale, T_shape>::value) {
    for (size_t n = 0; n < size(y_min); n++) {
      log_y_min[n] = log(value_of(y_min_vec[n]));
    }
  }

  VectorBuilder<include_summand<propto, T_shape>::value, T_partials_return,
                T_shape>
      log_alpha(size(alpha));
  if (include_summand<propto, T_shape>::value) {
    for (size_t n = 0; n < size(alpha); n++) {
      log_alpha[n] = log(value_of(alpha_vec[n]));
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    if (include_summand<propto, T_shape>::value) {
      logp += log_alpha[n];
    }
    if (include_summand<propto, T_scale, T_shape>::value) {
      logp += alpha_dbl * log_y_min[n];
    }
    if (include_summand<propto, T_y, T_shape>::value) {
      logp -= alpha_dbl * log_y[n] + log_y[n];
    }

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= alpha_dbl * inv_y[n] + inv_y[n];
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge2_.partials_[n] += alpha_dbl / value_of(y_min_vec[n]);
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge3_.partials_[n]
          += 1 / alpha_dbl + log_y_min[n] - log_y[n];
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_scale, T_shape> pareto_lpdf(const T_y& y,
                                                        const T_scale& y_min,
                                                        const T_shape& alpha) {
  return pareto_lpdf<false>(y, y_min, alpha);
}

}  // namespace math
}  // namespace stan
#endif
