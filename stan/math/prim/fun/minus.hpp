#ifndef STAN_MATH_PRIM_FUN_MINUS_HPP
#define STAN_MATH_PRIM_FUN_MINUS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>

namespace stan {
namespace math {

/**
 * Returns the negation of the specified scalar or matrix.
 *
 * @tparam T Type of subtrahend.
 * @param x Subtrahend.
 * @return Negation of subtrahend.
 */
template <typename T, require_not_std_vector_t<T>* = nullptr>
inline auto minus(T&& x) {
  return make_holder([](auto&& x_) { return -x_; }, std::forward<T>(x));
}

/**
 * Return the negation of the each element of a vector
 *
 * @tparam T Type of container.
 * @param x Container.
 * @return Container where each element is negated.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto minus(T&& x) {
  return apply_vector_unary<T>::apply(std::forward<T>(x),
                                      [](auto&& v) { return -v; });
}

}  // namespace math
}  // namespace stan

#endif
