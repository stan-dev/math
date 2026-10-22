#ifndef STAN_MATH_FWD_FUN_LOG_SOFTMAX_HPP
#define STAN_MATH_FWD_FUN_LOG_SOFTMAX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log_softmax.hpp>
#include <stan/math/prim/fun/softmax_log_sum_exp.hpp>
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
inline plain_type_t<Mat> log_softmax(Mat&& x) {
  if (x.size() == 0) {
    return {};
  }
  decltype(auto) x_ref = to_ref(std::forward<Mat>(x));
  const auto x_val = value_of(x_ref).eval();
  const auto d_in = x_ref.d().eval();
  const auto softmax_lse = internal::softmax_log_sum_exp(x_val);
  return to_fvar(x_val.array() - softmax_lse.second,
                 d_in.array() - softmax_lse.first.cwiseProduct(d_in).sum());
}

}  // namespace math
}  // namespace stan
#endif
