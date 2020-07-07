#ifndef STAN_MATH_PRIM_PROB_PARETO_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_PARETO_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_scale, typename T_shape>
return_type_t<T_y, T_scale, T_shape> pareto_lccdf(const T_y& y,
                                                  const T_scale& y_min,
                                                  const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_scale, T_shape>;
  using std::exp;
  using std::log;
  static const char* function = "pareto_lccdf";
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_positive_finite(function, "Scale parameter", y_min);
  check_positive_finite(function, "Shape parameter", alpha);
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         y_min, "Shape parameter", alpha);

  if (size_zero(y, y_min, alpha)) {
    return 0;
  }

  T_partials_return P(0.0);
  operands_and_partials<T_y, T_scale, T_shape> ops_partials(y, y_min, alpha);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_scale> y_min_vec(y_min);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  size_t N = max_size(y, y_min, alpha);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(y); i++) {
    if (value_of(y_vec[i]) < value_of(y_min_vec[i])) {
      return ops_partials.build(0.0);
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(y_vec[n]) == INFTY) {
      return ops_partials.build(negative_infinity());
    }

    const T_partials_return log_dbl
        = log(value_of(y_min_vec[n]) / value_of(y_vec[n]));
    const T_partials_return y_min_inv_dbl = 1.0 / value_of(y_min_vec[n]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);

    P += alpha_dbl * log_dbl;

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          -= alpha_dbl * y_min_inv_dbl * exp(log_dbl);
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge2_.partials_[n] += alpha_dbl * y_min_inv_dbl;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge3_.partials_[n] += log_dbl;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
