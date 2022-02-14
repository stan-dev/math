#ifndef STAN_MATH_REV_AGGREGATE_PARTIAL_HPP
#define STAN_MATH_REV_AGGREGATE_PARTIAL_HPP

#include <stan/math/rev/meta.hpp>

namespace stan {
namespace math {
namespace internal {
/*
template <typename T>
using base_matrix_t = std::conditional_t<
  is_var_matrix<std::decay_t<T>>::value,
  value_type_t<std::decay_t<T>>,
  std::decay_t<T>>;*/
} // namespace internal

template <
  typename ArgT, typename ReturnT, typename PartialT,
  require_eigen_col_vector_t<ArgT>* = nullptr,
  require_eigen_matrix_dynamic_t<ReturnT>* = nullptr,
  require_eigen_matrix_dynamic_t<PartialT>* = nullptr,
  require_st_var<return_type_t<ArgT, ReturnT, PartialT>>* = nullptr
>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, PartialT&& x) {
  return (rtn.adj().array() * x.array()).rowwise().sum().eval();
}

template <
  typename ArgT, typename ReturnT, typename PartialT,
  require_eigen_matrix_dynamic_t<ArgT>* = nullptr,
  require_eigen_matrix_dynamic_t<ReturnT>* = nullptr,
  require_eigen_col_vector_t<PartialT>* = nullptr,
  require_st_var<return_type_t<ArgT, ReturnT, PartialT>>* = nullptr
>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, PartialT&& x) {
  return (rtn.adj().array().colwise() * x.array()).eval();
}

}  // namespace math
}  // namespace stan
#endif
