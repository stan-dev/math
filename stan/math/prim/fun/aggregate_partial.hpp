#ifndef STAN_MATH_PRIM_FUN_AGGREGATE_PARTIAL_HPP
#define STAN_MATH_PRIM_FUN_AGGREGATE_PARTIAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/sum.hpp>

namespace stan {
namespace math {

template <typename ReturnT, typename PartialT,
          require_all_stan_scalar_t<ReturnT, PartialT>* = nullptr>
inline PartialT aggregate_partial(ReturnT&& rtn, const PartialT& x) {
  return x;
}

template <typename ReturnT, typename PartialT,
          require_all_container_or_var_matrix_t<ReturnT, PartialT>* = nullptr>
inline PartialT aggregate_partial(ReturnT&& rtn, const PartialT& x) {
  return x;
}

template <typename ReturnT, typename PartialT,
          require_stan_scalar_t<ReturnT>* = nullptr,
          require_not_stan_scalar_t<PartialT>* = nullptr>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, const PartialT& x) {
  return sum(x);
}

template <typename ReturnT, typename PartialT,
          require_eigen_t<ReturnT>* = nullptr,
          require_stan_scalar_t<PartialT>* = nullptr>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, const PartialT& x) {
  return plain_type_t<promote_scalar_t<PartialT, ReturnT>>::Constant(rtn.rows(), rtn.cols(), x).eval();
}


template <typename ReturnT, typename PartialT,
          require_all_var_matrix_t<ReturnT>* = nullptr,
          require_stan_scalar_t<PartialT>* = nullptr>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, const PartialT& x) {
  return plain_type_t<value_type_t<ReturnT>>::Constant(rtn.rows(), rtn.cols(), x).eval();
}

}  // namespace math
}  // namespace stan
#endif
