#ifndef STAN_MATH_PRIM_SCAL_PROB_STD_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_STD_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>

namespace stan {
namespace math {

/**
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
inline auto std_normal_lpdf(const T_y& y) {
  using T_partials = partials_return_t<T_y>;
  T_partials logp(0.0);

  static const char* function = "std_normal_lpdf";
  check_not_nan(function, "Random variable", y);

  const scalar_seq_view<T_y> y_vec(y);
  const auto size_y = length(y);
  operands_and_partials<T_y> ops_partials(y);
  if (!include_summand<propto, T_y>::value) {
    return ops_partials.build(logp);
  } else if (size_zero(y)) {
    return ops_partials.build(logp);
  }
  for (size_t n = 0; n < size_y; n++) {
    const T_partials y_val = value_of(y_vec[n]);
    logp += y_val * y_val;
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= y_val;
    }
  }
  logp *= -0.5;
  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * size_y;
  }
  return ops_partials.build(logp);
}

template <typename T_y>
inline auto std_normal_lpdf(const T_y& y) {
  return std_normal_lpdf<false>(y);
}

}  // namespace math
}  // namespace stan
#endif
