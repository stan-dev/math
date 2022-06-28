#ifndef STAN_MATH_PRIM_FUN_TO_INT_HPP
#define STAN_MATH_PRIM_FUN_TO_INT_HPP

#include <stan/math/prim/functor/apply_scalar_unary.hpp>

namespace stan {
namespace math {

/**
 * Returns the input scalar as an integer type
 *
 * @tparam T type of argument (must be arithmetic)
 * @param x argument
 * @return iIteger value of argument
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline int to_int(T x) {
  if (x < std::numeric_limits<int>::min() ||
      x > std::numeric_limits<int>::max()) {
    std::ostringstream msg;
    msg << "Value " << x << " is too large to be represented as an integer";
    throw std::invalid_argument(msg.str());
  }
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
