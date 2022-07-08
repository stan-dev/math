#ifndef STAN_MATH_PRIM_FUN_TO_INT_HPP
#define STAN_MATH_PRIM_FUN_TO_INT_HPP

#include <stan/math/prim/err/check_bounded.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>

namespace stan {
namespace math {

/**
 * Returns the input scalar as an integer type. Specialisation for integral
 * types which do not need conversion, reduces to a no-op.
 *
 * @tparam T type of integral argument
 * @param x argument
 * @return Input argument unchanged
 */
template <typename T, require_integral_t<T>* = nullptr>
inline T to_int(T x) {
  return std::forward<T>(x);
}

/**
 * Returns the input scalar as an integer type. This function performs no
 * rounding and simply truncates the decimal to return only the signficand as an
 * integer.
 *
 * Casting NaN and Inf values to integers is considered undefined behavior as
 * NaN and Inf cannot be represented as an integer and most implementations
 * simply overflow, as such this function throws for these inputs.
 *
 * The function also throws for floating-point values that are too large to be
 * represented as an integer.
 *
 * @tparam T type of argument (must be arithmetic)
 * @param x argument
 * @return Integer value of argument
 * @throw std::domain_error for NaN, Inf, or floating point values not in range
 *        to be represented as int
 */
template <typename T, require_floating_point_t<T>* = nullptr>
inline int to_int(T x) {
  static const char* function = "to_int";
  check_bounded(function, "x", x, std::numeric_limits<int>::min(),
                std::numeric_limits<int>::max());
  return static_cast<int>(x);
}

/**
 * Return elementwise integer value of the specified real-valued
 * container.
 *
 * @tparam T type of argument
 * @param x argument
 * @return Integer value of argument
 */
struct to_int_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return to_int(x);
  }
};

/**
 * Returns the elementwise `to_int()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x argument
 * @return Integer value of each variable in the container.
 */
template <typename Container,
          require_std_vector_st<std::is_arithmetic, Container>* = nullptr>
inline auto to_int(const Container& x) {
  return apply_scalar_unary<to_int_fun, Container>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
