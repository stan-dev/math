#ifndef STAN_MATH_PRIM_PROB_POISSON_LPMF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

// Poisson(n|lambda)  [lambda > 0;  n >= 0]
template <bool propto, typename T_n, typename T_rate,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_rate>* = nullptr>
return_type_t<T_rate> poisson_lpmf(const T_n& n, const T_rate& lambda) {
  using T_partials_return = partials_return_t<T_n, T_rate>;
  using T_n_ref = ref_type_if_t<!is_constant<T_n>::value, T_n>;
  using T_lambda_ref = ref_type_if_t<!is_constant<T_rate>::value, T_rate>;
  using std::isinf;
  static const char* function = "poisson_lpmf";
  check_consistent_sizes(function, "Random variable", n, "Rate parameter",
                         lambda);

  T_n_ref n_ref = n;
  T_lambda_ref lambda_ref = lambda;

  decltype(auto) n_val = to_ref(as_value_column_array_or_scalar(n_ref));
  decltype(auto) lambda_val
      = to_ref(as_value_column_array_or_scalar(lambda_ref));

  check_nonnegative(function, "Random variable", n_val);
  check_nonnegative(function, "Rate parameter", lambda_val);

  if (size_zero(n, lambda)) {
    return 0.0;
  }
  if (!include_summand<propto, T_rate>::value) {
    return 0.0;
  }
  if (sum(promote_scalar<int>(isinf(lambda_val)))) {
    return LOG_ZERO;
  }

  size_t N = max_size(n, lambda);
  scalar_seq_view<decltype(n_val)> n_vec(n_val);
  scalar_seq_view<decltype(lambda_val)> lambda_vec(lambda_val);
  for (size_t i = 0; i < N; i++) {
    if (lambda_vec[i] == 0 && n_vec[i] != 0) {
      return LOG_ZERO;
    }
  }

  operands_and_partials<T_lambda_ref> ops_partials(lambda_ref);

  T_partials_return logp = stan::math::sum(multiply_log(n_val, lambda_val));
  if (include_summand<propto, T_rate>::value) {
    logp -= sum(lambda_val) * N / math::size(lambda);
  }
  if (include_summand<propto>::value) {
    logp -= sum(lgamma(n_val + 1.0)) * N / math::size(n);
  }

  if (!is_constant_all<T_rate>::value) {
    ops_partials.edge1_.partials_ = n_val / lambda_val - 1.0;
  }

  return ops_partials.build(logp);
}

template <typename T_n, typename T_rate>
inline return_type_t<T_rate> poisson_lpmf(const T_n& n, const T_rate& lambda) {
  return poisson_lpmf<false>(n, lambda);
}

}  // namespace math
}  // namespace stan
#endif
