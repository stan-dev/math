#ifndef STAN_MATH_PRIM_PROB_POISSON_CDF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/gamma_q.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

// Poisson CDF
template <typename T_n, typename T_rate>
return_type_t<T_rate> poisson_cdf(const T_n& n, const T_rate& lambda) {
  static const char* function = "poisson_cdf";
  using T_partials_return = partials_return_t<T_n, T_rate>;

  if (size_zero(n, lambda)) {
    return 1.0;
  }

  T_partials_return P(1.0);

  check_not_nan(function, "Rate parameter", lambda);
  check_nonnegative(function, "Rate parameter", lambda);
  check_consistent_sizes(function, "Random variable", n, "Rate parameter",
                         lambda);

  using std::exp;
  using std::pow;

  operands_and_partials<T_rate> ops_partials(lambda);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_rate> lambda_vec(lambda);
  size_t size_n = size(n);
  size_t size_lambda = size(lambda);
  size_t max_size_seq_view = max_size(n, lambda);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size_n; i++) {
    if (value_of(n_vec[i]) < 0) {
      return ops_partials.build(0.0);
    }
  }

  VectorBuilder<!is_constant_all<T_rate>::value, T_partials_return, T_n>
      tgamma_n_plus_one(size_n);
  VectorBuilder<!is_constant_all<T_rate>::value, T_partials_return, T_rate>
      exp_minus_lambda(size_lambda);

  if (!is_constant_all<T_rate>::value) {
    for (size_t i = 0; i < size_n; i++) {
      tgamma_n_plus_one[i] = tgamma(n_vec[i] + 1.0);
    }
    for (size_t i = 0; i < size_lambda; i++) {
      exp_minus_lambda[i] = exp(-value_of(lambda_vec[i]));
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(n_vec[i]) == std::numeric_limits<int>::max()) {
      continue;
    }

    const T_partials_return n_dbl = value_of(n_vec[i]);
    const T_partials_return lambda_dbl = value_of(lambda_vec[i]);
    const T_partials_return Pi = gamma_q(n_dbl + 1, lambda_dbl);

    P *= Pi;

    if (!is_constant_all<T_rate>::value) {
      ops_partials.edge1_.partials_[i] -= exp_minus_lambda[i]
                                          * pow(lambda_dbl, n_dbl)
                                          / (tgamma_n_plus_one[i] * Pi);
    }
  }

  if (!is_constant_all<T_rate>::value) {
    for (size_t i = 0; i < size_lambda; ++i) {
      ops_partials.edge1_.partials_[i] *= P;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
