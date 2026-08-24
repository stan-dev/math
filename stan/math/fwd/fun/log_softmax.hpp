#ifndef STAN_MATH_FWD_FUN_LOG_SOFTMAX_HPP
#define STAN_MATH_FWD_FUN_LOG_SOFTMAX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun/softmax.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log_softmax.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>

namespace stan {
namespace math {

/**
 * Return the log softmax of each vector or matrix in a container of `fvar`
 * values.
 *
 * @tparam T `std::vector` whose scalar type is `fvar`
 * @param x container of vectors or matrices to transform
 * @return container of log softmax results
 */
template <typename T, require_std_vector_st<is_fvar, T>* = nullptr>
inline auto log_softmax(T&& x) {
  return apply_vector_unary<T>::apply(std::forward<T>(x), [](auto&& v) {
    return log_softmax(std::forward<decltype(v)>(v));
  });
}

/**
 * Return the log softmax of the specified vector or matrix of `fvar` values.
 *
 * @tparam Mat Eigen vector or matrix with `fvar` scalar
 * @param x vector or matrix to transform
 * @return log softmax of the vector or matrix, or an empty result if the
 * input is empty
 */
template <typename Mat, require_eigen_vt<is_fvar, Mat>* = nullptr>
inline auto log_softmax(Mat&& x) {
  using mat = std::decay_t<Mat>;
  constexpr int Rows = mat::RowsAtCompileTime;
  constexpr int Cols = mat::ColsAtCompileTime;
  using T = typename value_type_t<mat>::Scalar;
  decltype(auto) x_ref = to_ref(std::forward<Mat>(x));
  if (x_ref.size() == 0) {
    return Eigen::Matrix<fvar<T>, Rows, Cols>{};
  }
  const auto x_val = value_of(x_ref);
  const auto lse = log_sum_exp(x_val);
  const auto s = softmax(x_val);
  const auto d_in = x_ref.d();
  const auto dot_sd = (s.array() * d_in.array()).sum();

  Eigen::Matrix<fvar<T>, Rows, Cols> result(x_ref.rows(), x_ref.cols());
  result.val() = (x_val.array() - lse).matrix();
  result.d() = (d_in.array() - dot_sd).matrix();
  return result;
}

}  // namespace math
}  // namespace stan
#endif