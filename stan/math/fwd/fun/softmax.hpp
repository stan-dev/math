#ifndef STAN_MATH_FWD_FUN_SOFTMAX_HPP
#define STAN_MATH_FWD_FUN_SOFTMAX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
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
 * Return the softmax of the specified vector of `fvar` values.
 *
 * @tparam Mat Eigen vector or matrix with `fvar` scalar
 * @param x vector or matrix to transform
 * @return softmax of the vector, matrix, or an empty result if the input is
 *   empty
 */
template <typename Mat, require_eigen_vt<is_fvar, Mat>* = nullptr>
inline plain_type_t<Mat> softmax(Mat&& x) {
  if (x.size() == 0) {
    return {};
  }
  decltype(auto) x_ref = to_ref(std::forward<Mat>(x));
  const auto s = softmax(value_of(x_ref)).eval();
  const auto d_in = x_ref.d().eval();

  return to_fvar(s, s.array() * (d_in.array() - s.cwiseProduct(d_in).sum()));
}

}  // namespace math
}  // namespace stan
#endif
