#ifndef STAN_MATH_REV_CORE_OPERATOR_EQUAL_HPP
#define STAN_MATH_REV_CORE_OPERATOR_EQUAL_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Equality operator comparing two variables' values (C++).
 *
   \f[
   \mbox{operator==}(x, y) =
   \begin{cases}
     0 & \mbox{if } x \neq y\\
     1 & \mbox{if } x = y \\[6pt]
     0 & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a First variable.
 * @param b Second variable.
 * @return True if the first variable's value is the same as the
 * second's.
 */
inline bool operator==(const var& a, const var& b) {
  return a.val() == b.val();
}

/**
 * Equality operator comparing a variable's value and a double
 * (C++).
 *
 * @tparam Arith An arithmetic type
 * @param a First variable.
 * @param b Second value.
 * @return True if the first variable's value is the same as the
 * second value.
 */
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline bool operator==(const var& a, Arith b) {
  return a.val() == b;
}

/**
 * Equality operator comparing a scalar and a variable's value
 * (C++).
 *
 * @tparam Arith An arithmetic type
 * @param a First scalar.
 * @param b Second variable.
 * @return True if the variable's value is equal to the scalar.
 */
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline bool operator==(Arith a, const var& b) {
  return a == b.val();
}

/**
 * Return `true` if the real number is equal to the real part of the
 * complex number, and the imaginary part of the complex number is zero
 *
 * @param x real number
 * @param z complex number
 * @return `true` if the real number is equal to the real part of the
 * complex number, and the imaginary part of the complex number is zero
 */
inline bool operator==(const var& x, const std::complex<var>& z) {
  return x == z.real() && 0 == z.imag();
}

/**
 * Return `true` if the real number is equal to the real part of the
 * complex number, and the imaginary part of the complex number is zero
 *
 * @param z complex number
 * @param y real number
 * @return `true` if the real number is equal to the real part of the
 * complex number, and the imaginary part of the complex number is zero
 */
inline bool operator==(const std::complex<var>& z, const var& y) {
  return z.real() == y && z.imag() == 0;
}

}  // namespace math
}  // namespace stan
#endif
