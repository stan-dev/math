#ifndef STAN_MATH_FWD_FUN_AGGREGATE_PARTIAL_HPP
#define STAN_MATH_FWD_FUN_AGGREGATE_PARTIAL_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun/Eigen_NumTraits.hpp>

namespace stan {
namespace math {

template <
  typename ArgT, typename ReturnT, typename PartialT,
  require_eigen_matrix_dynamic_t<ArgT>* = nullptr,
  require_eigen_col_vector_t<ReturnT>* = nullptr,
  require_eigen_matrix_dynamic_t<PartialT>* = nullptr,
  require_st_fvar<return_type_t<ArgT, ReturnT, PartialT>>* = nullptr
>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, PartialT&& x) {
  return elt_multiply(rtn.d().replicate(1, x.cols()), x).array().eval();
}

template <
  typename ArgT, typename ReturnT, typename PartialT,
  require_eigen_matrix_dynamic_t<ArgT>* = nullptr,
  require_eigen_matrix_dynamic_t<ReturnT>* = nullptr,
  require_eigen_col_vector_t<PartialT>* = nullptr,
  require_st_fvar<return_type_t<ArgT, ReturnT, PartialT>>* = nullptr
>
inline decltype(auto) aggregate_partial(ReturnT&& rtn, PartialT&& x) {
  return elt_multiply(rtn.d(), x.replicate(1, rtn.cols())).array().eval();
}

}  // namespace math
}  // namespace stan
#endif
