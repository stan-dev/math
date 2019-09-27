#ifndef STAN_MATH_PRIM_SCAL_PROB_POISSON_LCCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_POISSON_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/gamma_p.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

template <typename T_n, typename T_rate>
return_type_t<T_rate> poisson_lccdf(const T_n& n, const T_rate& lambda) {
  static const char* function = "poisson_lccdf";
  using T_partials_return = partials_return_t<T_n, T_rate>;

  if (size_zero(n, lambda)) {
    return 0.0;
  }

  T_partials_return P(0.0);

  check_not_nan(function, "Rate parameter", lambda);
  check_nonnegative(function, "Rate parameter", lambda);
  check_consistent_sizes(function, "Random variable", n, "Rate parameter",
                         lambda);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_rate> lambda_vec(lambda);
  size_t size = max_size(n, lambda);

  using std::exp;
  using std::log;

  operands_and_partials<T_rate> ops_partials(lambda);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as neg infinity
  for (size_t i = 0; i < stan::length(n); i++) {
    if (value_of(n_vec[i]) < 0) {
      return ops_partials.build(0.0);
    }
  }

  for (size_t i = 0; i < size; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(n_vec[i]) == std::numeric_limits<int>::max()) {
      return ops_partials.build(negative_infinity());
    }

    const T_partials_return n_dbl = value_of(n_vec[i]);
    const T_partials_return lambda_dbl = value_of(lambda_vec[i]);
    const T_partials_return log_Pi = log(gamma_p(n_dbl + 1, lambda_dbl));

    P += log_Pi;

    if (!is_constant_all<T_rate>::value) {
      ops_partials.edge1_.partials_[i] += exp(
          n_dbl * log(lambda_dbl) - lambda_dbl - lgamma(n_dbl + 1) - log_Pi);
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
