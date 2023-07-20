#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the normal density for the specified scalar(s) given
 * a location of 0 and a scale of 1. y can be either
 * a scalar or a vector.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation.
 *
 * @tparam T_y type of scalar
 * @param y (Sequence of) scalar(s).
 * @return The log of the product of the densities.
 * @throw std::domain_error if any scalar is nan.
 */
template <
    bool propto, typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>* = nullptr>
return_type_t<T_y> std_normal_lpdf(const T_y& y) {
  using T_partials_return = partials_return_t<T_y>;
  using T_y_ref = ref_type_t<T_y>;
  static constexpr const char* function = "std_normal_lpdf";
  T_y_ref y_ref = y;
  check_not_nan(function, "Random variable", y_ref);

  if (size_zero(y)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y>::value) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  auto ops_partials = make_partials_propagator(y_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  size_t N = stan::math::size(y);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_val = y_vec.val(n);
    logp += y_val * y_val;
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials)[n] -= y_val;
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
