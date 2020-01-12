#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the normal density for the specified scalar(s) given
 * a location of 0 and a scale of 1. y can be either
 * a scalar or a vector.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation.
 * @tparam T_y Underlying type of scalar in sequence.
 * @param y (Sequence of) scalar(s).
 * @return The log of the product of the densities.
 * @throw std::domain_error if any scalar is nan.
 */
template <bool propto, typename T_y>
return_type_t<T_y> std_normal_lpdf(const T_y& y) {
  static const char* function = "std_normal_lpdf";
  using T_partials_return = partials_return_t<T_y>;

  if (size_zero(y)) {
    return 0.0;
  }

  check_not_nan(function, "Random variable", y);

  if (!include_summand<propto, T_y>::value) {
    return 0.0;
  }

  operands_and_partials<T_y> ops_partials(y);

  T_partials_return logp(0.0);
  scalar_seq_view<T_y> y_vec(y);
  size_t N = size(y);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_val = value_of(y_vec[n]);
    logp += y_val * y_val;
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= y_val;
    }
  }
  logp *= -0.5;
  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * N;
  }

  return ops_partials.build(logp);
}

template <typename T_y>
inline return_type_t<T_y> std_normal_lpdf(const T_y& y) {
  return std_normal_lpdf<false>(y);
}

}  // namespace math
}  // namespace stan
#endif
