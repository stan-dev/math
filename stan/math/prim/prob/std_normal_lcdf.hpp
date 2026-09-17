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
#include <stan/math/prim/prob/std_normal_lcdf_impl.hpp>

namespace stan {
namespace math {
namespace internal {
constexpr char std_normal_lcdf_func[] = "std_normal_lcdf";
}  // namespace internal

/** \ingroup prob_dists
 * @brief Calculates the log of the cdf of the standard normal distribution
 *
 * Shares the scalar value and slope calculation with normal_lcdf through
 * std_normal_lcdf_impl.hpp.
 *
 * @tparam func name reported by the error checks. Reflected distributions
 *   such as `std_normal_lccdf` delegate here and pass their own name so that
 *   exceptions name the function the user actually called.
 * @tparam T_y A vector or scalar type for the random variable.
 * @param y (Sequence of) scalar(s).
 * @return The log of the standard normal cdf evaluated at the specified
 *   argument. If given a container, the log of the product of the cdfs.
 */
template <
    const char* func = internal::std_normal_lcdf_func, typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>* = nullptr>
inline return_type_t<T_y> std_normal_lcdf(const T_y& y) {
  using T_y_ref = ref_type_t<T_y>;
  static constexpr const char* function = func;
  const T_y_ref y_ref = y;
  check_not_nan(function, "Random variable", y_ref);

  if (size_zero(y_ref)) {
    return 0;
  }

  const auto& y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  const auto [values, slopes]
      = internal::std_normal_lcdf_value_grad<is_autodiff_v<T_y>>(y_val);
  auto ops_partials = make_partials_propagator(y_ref);
  const auto log_p = sum(values);
  if constexpr (is_autodiff_v<T_y>) {
    partials<0>(ops_partials) = slopes;
  }
  return ops_partials.build(log_p);
}

}  // namespace math
}  // namespace stan
#endif
