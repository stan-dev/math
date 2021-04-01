#ifndef STAN_MATH_REV_FUN_ADJOINT_OF_HPP
#define STAN_MATH_REV_FUN_ADJOINT_OF_HPP

#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

namespace internal {
struct nonexisting_adjoint {
  template <typename T>
  nonexisting_adjoint operator+(const T&) {
    return *this;
  }
  template <typename T>
  nonexisting_adjoint operator+=(T) {
    throw std::runtime_error(
        "internal::nonexisting_adjoint::operator+= should never be called! "
        "Please file a bug report.");
  }
  template <typename T>
  nonexisting_adjoint operator-=(T) {
    throw std::runtime_error(
        "internal::nonexisting_adjoint::operator-= should never be called! "
        "Please file a bug report.");
  }
};
}  // namespace internal

/**
 * Returns a reference to a variable's adjoint.
 *
 * @param x a var
 * @return reference to `x`'s adjoint
 */
template <typename T, require_var_t<T>* = nullptr>
auto& adjoint_of(const T& x) {
  return x.adj();
}

/**
 * Returns a reference to a variable's adjoint. If the input object is not var,
 * it does not have an adjoint and this returns a dummy object. It defines
 * operators += and -=, but they should not actually be called.
 *
 * @param x any non-var object
 * @return a dummy adjoint
 */
template <typename T, require_not_var_t<T>* = nullptr>
internal::nonexisting_adjoint adjoint_of(const T& x) {
  return {};
}

}  // namespace math
}  // namespace stan

#endif  // ADJOINT_OF_HPP
