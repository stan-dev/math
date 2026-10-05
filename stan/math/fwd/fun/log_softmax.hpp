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
 * Return the log softmax of the specified vector of `fvar` values.
 *
 * @tparam Vec Eigen vector with `fvar` scalar
 * @param x vector to transform
 * @return log softmax of the vector, or an empty result if the input is empty
 */
template <typename Vec, require_eigen_vector_vt<is_fvar, Vec>* = nullptr>
inline auto log_softmax(Vec&& x) {
  using vec = std::decay_t<Vec>;
  constexpr int Rows = vec::RowsAtCompileTime;
  constexpr int Cols = vec::ColsAtCompileTime;
  using T = typename value_type_t<vec>::Scalar;
  decltype(auto) x_ref = to_ref(std::forward<Vec>(x));
  if (x_ref.size() == 0) {
    return Eigen::Matrix<fvar<T>, Rows, Cols>{};
  }
  const auto s = softmax(value_of(x_ref));
  const auto d_in = x_ref.d();
  const auto dot_sd = s.dot(d_in);
  Eigen::Matrix<fvar<T>, Rows, Cols> result(x_ref.size());
  result.val() = s.array().log().matrix();
  result.d() = (d_in.array() - dot_sd).matrix();
  return result;
}

}  // namespace math
}  // namespace stan
#endif
