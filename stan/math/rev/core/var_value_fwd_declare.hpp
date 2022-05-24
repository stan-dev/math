#ifndef STAN_MATH_REV_CORE_VAR_VALUE_FWD_DECLARE_HPP
#define STAN_MATH_REV_CORE_VAR_VALUE_FWD_DECLARE_HPP

namespace stan {
namespace math {
/**
 * Independent (input) and dependent (output) variables for gradients.
 *
 * This class acts as a smart pointer, with resources managed by
 * an arena-based memory manager scoped to a single gradient
 * calculation.
 *
 * A var is constructed with a type `T` and used like any
 * other scalar. Arithmetical functions like negation, addition,
 * and subtraction, as well as a range of mathematical functions
 * like exponentiation and powers are overridden to operate on
 * var values objects.
 * @tparam T An Floating point type.
 */
template <typename T, typename = void>
class var_value;
}  // namespace math
}  // namespace stan
#endif
