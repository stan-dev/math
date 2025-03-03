#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/as_value_column_vector_or_scalar.hpp>
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

  const auto& y_val = as_value_column_vector_or_scalar(y_ref);
  T_partials_return logp = -dot_self(y_val) / 2.0;
  auto ops_partials = make_partials_propagator(y_ref);

  if (!is_constant_all<T_y>::value) {
    partials<0>(ops_partials) = -y_val;
  }

  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * math::size(y);
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
