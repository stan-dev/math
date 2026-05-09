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
 * Return the log softmax of the specified row vector of `fvar` values.
 * Delegates to the column-vector overload via transposition.
 *
 * @tparam RowVec Eigen row vector with `fvar` scalar
 * @param x row vector to transform
 * @return log softmax of the row vector
 */
template <typename RowVec, require_eigen_row_vector_t<RowVec>* = nullptr,
          require_t<is_fvar<value_type_t<RowVec>>>* = nullptr>
inline auto log_softmax(const RowVec& x) {
  return log_softmax(x.transpose()).transpose().eval();
}

/**
 * Return the log softmax of each vector in a container of `fvar` values.
 *
 * @tparam T `std::vector` whose scalar type is `fvar`
 * @param x container of vectors to transform
 * @return container of log softmax results
 */
template <typename T, require_std_vector_st<is_fvar, T>* = nullptr>
inline auto log_softmax(T&& x) {
  return apply_vector_unary<T>::apply(std::forward<T>(x), [](auto&& v) {
    return log_softmax(std::forward<decltype(v)>(v));
  });
}

/**
 * Return the log softmax of the specified column vector of `fvar` values.
 *
 * @tparam ColVec Eigen column vector with `fvar` scalar
 * @param x column vector to transform
 * @return log softmax of the column vector
 * @throw std::domain_error if the input size is 0
 */
template <typename ColVec,
          require_eigen_col_vector_vt<is_fvar, ColVec>* = nullptr>
inline auto log_softmax(const ColVec& x) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using T = typename value_type_t<ColVec>::Scalar;
  check_nonzero_size("log_softmax", "x", x);
  const auto& x_ref = to_ref(x);
  const Matrix<T, Dynamic, 1> s = softmax(value_of(x_ref));
  const auto d_in = x_ref.d().eval();
  const T dot_sd = s.dot(d_in);
  Matrix<fvar<T>, Dynamic, 1> result(x_ref.size());
  result.val() = s.array().log().matrix();
  result.d() = (d_in.array() - dot_sd).matrix();
  return result;
}

}  // namespace math
}  // namespace stan
#endif
