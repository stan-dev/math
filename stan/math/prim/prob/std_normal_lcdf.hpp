#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/fun/std_normal_lcdf_impl.hpp>

namespace stan {
namespace math {
namespace internal {

/** Log of the standard normal cdf, or of its complement when `reflect` is
 * set: the input is negated so no negated autodiff operands are built.
 */
template <bool reflect, typename T_y>
inline return_type_t<T_y> std_normal_lcdf_impl(const char* function, T_y&& y) {
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  constexpr double sign = reflect ? -1.0 : 1.0;
  T_y_ref y_ref = std::forward<T_y>(y);
  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  check_not_nan(function, "Random variable", y_val);

  if (size_zero(y_ref)) {
    return 0;
  }

  // Branching avoids materialising sign * y_val when not reflecting.
  const auto [values, slopes] = [&]() {
    if constexpr (reflect) {
      return std_normal_lcdf_value_grad<is_autodiff_v<T_y>>(-y_val);
    } else {
      return std_normal_lcdf_value_grad<is_autodiff_v<T_y>>(y_val);
    }
  }();
  auto ops_partials = make_partials_propagator(y_ref);
  if constexpr (is_autodiff_v<T_y>) {
    partials<0>(ops_partials) = sign * slopes;
  }
  return ops_partials.build(sum(values));
}

}  // namespace internal

/** \ingroup prob_dists
 * @brief Calculates the log of the cdf of the standard normal distribution
 *
 * Shares the scalar value and slope calculation with normal_lcdf through
 * std_normal_lcdf_impl.hpp.
 *
 * @tparam T_y A vector or scalar type for the random variable.
 * @param y (Sequence of) scalar(s).
 * @return The log of the standard normal cdf evaluated at the specified
 *   argument. If given a container, the log of the product of the cdfs.
 */
template <
    typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>* = nullptr>
inline return_type_t<T_y> std_normal_lcdf(T_y&& y) {
  return internal::std_normal_lcdf_impl<false>("std_normal_lcdf",
                                               std::forward<T_y>(y));
}

}  // namespace math
}  // namespace stan
#endif
