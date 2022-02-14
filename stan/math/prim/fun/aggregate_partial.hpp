#ifndef STAN_MATH_PRIM_FUN_AGGREGATE_PARTIAL_HPP
#define STAN_MATH_PRIM_FUN_AGGREGATE_PARTIAL_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {
namespace internal {

template <typename T>
using base_matrix_t = std::conditional_t<
  is_var_matrix<std::decay_t<T>>::value,
  value_type_t<std::decay_t<T>>,
  std::decay_t<T>>;
} // namespace internal

template <typename ArgT, typename ReturnT, typename PartialT,
          require_all_stan_scalar_t<ReturnT, PartialT>* = nullptr>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, PartialT&& x) {
  return rtn * x;
}

}  // namespace math
}  // namespace stan
#endif
