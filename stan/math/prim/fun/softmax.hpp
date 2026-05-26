#ifndef STAN_MATH_PRIM_FUN_SOFTMAX_HPP
#define STAN_MATH_PRIM_FUN_SOFTMAX_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the softmax of the specified vector.
 *
 * <p>
 * \f$
 * \mbox{softmax}(y)
 * = \frac{\exp(y)}
 * {\sum_{k=1}^K \exp(y_k)},
 * \f$
 *
 * <p>The entries in the Jacobian of the softmax function are given by
 * \f$
 * \begin{array}{l}
 * \displaystyle
 * \frac{\partial}{\partial y_m} \mbox{softmax}(y)[k]
 * \\[8pt]
 * \displaystyle
 * \mbox{ } \ \ \ = \left\{
 * \begin{array}{ll}
 * \mbox{softmax}(y)[k] \times (1 - \mbox{softmax}(y)[m])
 * & \mbox{ if } m = k, \mbox{ and}
 * \\[6pt]
 * -\mbox{softmax}(y)[k] \times \mbox{softmax}(y)[m]
 * & \mbox{ if } m \neq k.
 * \end{array}
 * \right.
 * \end{array}
 * \f$
 *
 * @tparam Vec type of the input vector
 * @param[in] v Vector to transform.
 * @return Unit simplex result of the softmax transform of the vector.
 */
template <typename Vec,
          require_eigen_vector_vt<std::is_arithmetic, Vec>* = nullptr>
inline plain_type_t<Vec> softmax(Vec&& v) {
  if (v.size() == 0) {
    return v;
  }
  decltype(auto) v_ref = to_ref(std::forward<Vec>(v));
  const auto theta = (v_ref.array() - v_ref.maxCoeff()).exp();
  return (theta / theta.sum()).matrix();
}

/**
 * Return the softmax of each vector in an array.
 *
 * @tparam T `std::vector` whose scalar type is arithmetic
 * @param[in] x Array of vectors to transform.
 * @return Array of unit simplex results.
 */
template <typename T, require_std_vector_st<std::is_arithmetic, T>* = nullptr>
inline auto softmax(T&& x) {
  return apply_vector_unary<T>::apply(std::forward<T>(x), [](auto&& v) {
    return softmax(std::forward<decltype(v)>(v));
  });
}

}  // namespace math
}  // namespace stan

#endif
