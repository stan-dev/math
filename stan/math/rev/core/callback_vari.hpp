#ifndef STAN_MATH_REV_CORE_CALLBACK_VARI_HPP
#define STAN_MATH_REV_CORE_CALLBACK_VARI_HPP

#include <stan/math/rev/core/vari.hpp>
#include <utility>

namespace stan {
namespace math {
namespace internal {

template <typename T, typename F>
struct callback_vari : public vari_value<T> {
  F rev_functor_;

  template <typename S,
            require_same_t<plain_type_t<T>, plain_type_t<S>>* = nullptr>
  explicit callback_vari(S&& value, F&& rev_functor)
      : vari_value<T>(std::move(value), true),
        rev_functor_(std::forward<F>(rev_functor)) {}

  inline void chain() final { rev_functor_(*this); }
};

}  // namespace internal

/**
 * Creates a new vari with given value and a callback that implements the
 * reverse pass (chain). The callback needs to accept a referenct to the vari.
 * If it needs any other data it should be implemented as a lambda capturing the
 * variables it needs.
 *
 * All captured values must be trivially destructible or they will leak memory!
 * `to_arena()` function can be used to ensure that.
 *
 * @tparam T type of value
 * @tparam F type of callable
 * @param value value of the vari
 * @param functor functor or other callable to call in the reverse pass
 */
template <typename T, typename F>
internal::callback_vari<plain_type_t<T>, F>* make_callback_vari(T&& value,
                                                                F&& functor) {
  return new internal::callback_vari<plain_type_t<T>, F>(
      std::move(value), std::forward<F>(functor));
}

/**
 * Creates a new var initialized with a callback_vari with a
 * given value and reverse-pass callback functor. The callback functor
 * will be passed a reference to the constructed vari.
 *
 * All captured values must be trivially destructible or they will leak memory.
 * `to_arena()` function can be used to ensure that.
 *
 * @tparam T type of value
 * @tparam F type of callable
 * @param value value of the vari
 * @param functor functor or other callable to call in the reverse pass
 */
template <typename T, typename F>
var_value<plain_type_t<T>> make_callback_var(T&& value, F&& functor) {
  return var_value<plain_type_t<T>>(
      make_callback_vari(std::move(value), std::forward<F>(functor)));
}

}  // namespace math
}  // namespace stan

#endif
