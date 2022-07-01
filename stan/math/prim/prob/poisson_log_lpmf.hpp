#ifndef STAN_MATH_PRIM_PROB_POISSON_LOG_LPMF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_LOG_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

// PoissonLog(n|alpha)  [n >= 0]   = Poisson(n|exp(alpha))
template <bool propto, typename T_n, typename T_log_rate,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_log_rate>* = nullptr>
return_type_t<T_log_rate> poisson_log_lpmf(const T_n& n,
                                           const T_log_rate& alpha) {
  using T_partials_return = partials_return_t<T_n, T_log_rate>;
  using T_n_ref = ref_type_if_t<!is_constant<T_n>::value, T_n>;
  using T_alpha_ref
      = ref_type_if_t<!is_constant<T_log_rate>::value, T_log_rate>;
  using std::isinf;
  static const char* function = "poisson_log_lpmf";
  check_consistent_sizes(function, "Random variable", n, "Log rate parameter",
                         alpha);

  T_n_ref n_ref = n;
  T_alpha_ref alpha_ref = alpha;

  decltype(auto) n_val = to_ref(as_value_column_array_or_scalar(n_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));

  check_nonnegative(function, "Random variable", n_val);
  check_not_nan(function, "Log rate parameter", alpha_val);

  if (size_zero(n, alpha)) {
    return 0.0;
  }
  if (!include_summand<propto, T_log_rate>::value) {
    return 0.0;
  }

  if (sum(promote_scalar<int>(INFTY == alpha_val))) {
    return LOG_ZERO;
  }

  size_t N = max_size(n, alpha);
  scalar_seq_view<decltype(n_val)> n_vec(n_val);
  scalar_seq_view<decltype(alpha_val)> alpha_vec(alpha_val);
  for (size_t i = 0; i < N; i++) {
    if (NEGATIVE_INFTY == alpha_vec[i] && n_vec[i] != 0) {
      return LOG_ZERO;
    }
  }

  operands_and_partials<T_alpha_ref> ops_partials(alpha_ref);

  const auto& exp_alpha
      = to_ref_if<!is_constant_all<T_log_rate>::value>(exp(alpha_val));

  T_partials_return logp = sum(n_val * alpha_val);
  if (include_summand<propto, T_log_rate>::value) {
    logp -= sum(exp_alpha) * N / math::size(alpha);
  }
  if (include_summand<propto>::value) {
    logp -= sum(lgamma(n_val + 1.0)) * N / math::size(n);
  }

  if (!is_constant_all<T_log_rate>::value) {
    ops_partials.edge1_.partials_ = n_val - exp_alpha;
  }

  return ops_partials.build(logp);
}

template <typename T_n, typename T_log_rate>
inline return_type_t<T_log_rate> poisson_log_lpmf(const T_n& n,
                                                  const T_log_rate& alpha) {
  return poisson_log_lpmf<false>(n, alpha);
}

}  // namespace math
}  // namespace stan
#endif
