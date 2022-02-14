#ifndef STAN_MATH_PRIM_FUN_AGGREGATE_PARTIAL_HPP
#define STAN_MATH_PRIM_FUN_AGGREGATE_PARTIAL_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename ArgT, typename ReturnT, typename PartialT,
          require_all_stan_scalar_t<ReturnT, PartialT>* = nullptr>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, PartialT&& x) {
  return rtn.adj() * x;
}

}  // namespace math
}  // namespace stan
#endif
