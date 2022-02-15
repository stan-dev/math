#ifndef STAN_MATH_REV_FUN_AGGREGATE_PARTIAL_HPP
#define STAN_MATH_REV_FUN_AGGREGATE_PARTIAL_HPP

#include <stan/math/rev/meta.hpp>

namespace stan {
namespace math {
namespace internal {
template <typename T>
using matrix_t = std::conditional_t<is_any_var_matrix<T>::value,
                                    value_type_t<std::decay_t<T>>,
                                    plain_type_t<std::decay_t<T>>>;
}

template <typename ArgT, typename ReturnT, typename PartialT,
          require_all_stan_scalar_t<ReturnT, PartialT>* = nullptr,
          require_st_var<return_type_t<ArgT, ReturnT, PartialT>>* = nullptr>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, PartialT&& x) {
  return rtn.adj() * x;
}
template <typename ArgT, typename ReturnT, typename PartialT,
          require_all_eigen_matrix_dynamic_t<
              internal::matrix_t<ArgT>, internal::matrix_t<ReturnT>,
              internal::matrix_t<PartialT>>* = nullptr,
          require_st_var<return_type_t<ArgT, ReturnT, PartialT>>* = nullptr>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, PartialT&& x) {
  return elt_multiply(rtn.adj(), x).array();
}
template <typename ArgT, typename ReturnT, typename PartialT,
          require_all_eigen_col_vector_t<
              internal::matrix_t<ArgT>, internal::matrix_t<ReturnT>,
              internal::matrix_t<PartialT>>* = nullptr,
          require_st_var<return_type_t<ArgT, ReturnT, PartialT>>* = nullptr>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, PartialT&& x) {
  return elt_multiply(rtn.adj(), x).array();
}
template <typename ArgT, typename ReturnT, typename PartialT,
          require_all_eigen_row_vector_t<
              internal::matrix_t<ArgT>, internal::matrix_t<ReturnT>,
              internal::matrix_t<PartialT>>* = nullptr,
          require_st_var<return_type_t<ArgT, ReturnT, PartialT>>* = nullptr>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, PartialT&& x) {
  return elt_multiply(rtn.adj(), x).array();
}

template <
  typename ArgT, typename ReturnT, typename PartialT,
  require_eigen_col_vector_t<internal::matrix_t<ArgT>>* = nullptr,
  require_eigen_matrix_dynamic_t<internal::matrix_t<ReturnT>>* = nullptr,
  require_eigen_matrix_dynamic_t<internal::matrix_t<PartialT>>* = nullptr,
  require_st_var<return_type_t<ArgT, ReturnT, PartialT>>* = nullptr
>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, PartialT&& x) {
  return (rtn.adj().array() * x.array()).rowwise().sum().eval();
}

template <
  typename ArgT, typename ReturnT, typename PartialT,
  require_eigen_matrix_dynamic_t<internal::matrix_t<ArgT>>* = nullptr,
  require_eigen_matrix_dynamic_t<internal::matrix_t<ReturnT>>* = nullptr,
  require_eigen_col_vector_t<internal::matrix_t<PartialT>>* = nullptr,
  require_st_var<return_type_t<ArgT, ReturnT, PartialT>>* = nullptr
>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, PartialT&& x) {
  return (rtn.adj().array().colwise() * x.array()).eval();
}

}  // namespace math
}  // namespace stan
#endif
