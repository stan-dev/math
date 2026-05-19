#ifndef STAN_MATH_FWD_FUN_SOFTMAX_HPP
#define STAN_MATH_FWD_FUN_SOFTMAX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/value_of.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/softmax.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>

namespace stan {
namespace math {

/**
 * Return the softmax of each vector in a container of `fvar` values.
 *
 * @tparam T `std::vector` whose scalar type is `fvar`
 * @param x container of vectors to transform
 * @return container of softmax results
 */
template <typename T, require_std_vector_st<is_fvar, T>* = nullptr>
inline auto softmax(T&& x) {
  return apply_vector_unary<T>::apply(std::forward<T>(x), [](auto&& v) {
    return softmax(std::forward<decltype(v)>(v));
  });
}

/**
 * Return the softmax of the specified column vector of `fvar` values.
 *
 * @tparam ColVec Eigen column vector with `fvar` scalar
 * @param x column vector to transform
 * @return softmax of the column vector
 */
template <typename Vec, require_eigen_vector_vt<is_fvar, Vec>* = nullptr>
inline auto softmax(Vec&& x) {
  using vec = std::decay_t<Vec>;
  constexpr int Rows = vec::RowsAtCompileTime;
  constexpr int Cols = vec::ColsAtCompileTime;
  using T = typename value_type_t<vec>::Scalar;
  if (x.size() == 0) {
    return Eigen::Matrix<fvar<T>, Rows, Cols>();
  }
  decltype(auto) x_ref = to_ref(std::forward<Vec>(x));
  const auto s = softmax(value_of(x_ref));
  const auto d_in = x_ref.d();
  const auto dot_sd = s.dot(d_in);
  Eigen::Matrix<fvar<T>, Rows, Cols> result(x_ref.size());
  result.val() = s;
  result.d() = (s.array() * (d_in.array() - dot_sd)).matrix();
  return result;
}

}  // namespace math
}  // namespace stan
#endif
