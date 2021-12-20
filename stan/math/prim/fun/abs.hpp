#ifndef STAN_MATH_PRIM_FUN_ABS_HPP
#define STAN_MATH_PRIM_FUN_ABS_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/hypot.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the absolute value of the specified arithmetic argument.
 * The return type is the same as the argument type.
 *
 * @tparam T type of arument (must be arithmetic)
 * @param x argument
 * @return absolute value of argument
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
T abs(T x) {
  return std::abs(x);
}

/*
 * Return the elementwise absolute value of the specified container.
 *
 * @tparam T type of elements in the vector
 * @param x vector argument
 * @return elementwise absolute value of argument
 */
template <typename T>
std::vector<T> abs(const std::vector<T>& x) {
  std::vector<T> y(x.size());
  for (size_t n = 0; n < x.size(); ++n)
    y[n] = abs(x[n]);
  return y;
}

/**
 * Return the elementwise absolute value of the specified matrix,
 * vector, or row vector.
 *
 * @tparam T type of scalar for matrix argument (real or complex)
 * @tparam R row specification (1 or -1)
 * @tparam C column specification (1 or -1)
 * @param x argument
 * @return elementwise absolute value of argument
 */
template <typename T, int R, int C>
Eigen::Matrix<T, R, C> abs(const Eigen::Matrix<T, R, C>& x) {
  return fabs(x);
}

/**
 * Return the absolute value (also known as the norm, modulus, or
 * magnitude) of the specified complex argument.
 *
 * @tparam T type of argument (must be complex)
 * @param x argument
 * @return absolute value of argument (a real number)
 */
template <typename T, require_complex_t<T>* = nullptr>
auto abs(T x) {
  return hypot(x.real(), x.imag());
}

/**
 * Return elementwise absolute value of the specified real-valued
 * container.
 *
 * @tparam T type of argument
 * @param x argument
 * @return absolute value of argument
 */
struct abs_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return fabs(x);
  }
};

/**
 * Returns the elementwise `abs()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x argument
 * @return Absolute value of each variable in the container.
 */
// template <typename Container,
//           require_not_container_st<std::is_arithmetic, Container>* = nullptr,
//           require_not_var_matrix_t<Container>* = nullptr,
//           require_not_stan_scalar_t<Container>* = nullptr>
// inline auto abs(const Container& x) {
//   return apply_scalar_unary<abs_fun, Container>::apply(x);
// }

// /**
//  * Version of `abs()` that accepts std::vectors, Eigen Matrix/Array objects
//  * or expressions, and containers of these.
//  *
//  * @tparam Container Type of x
//  * @param x argument
//  * @return Absolute value of each variable in the container.
//  */
// template <typename Container,
//           require_container_st<std::is_arithmetic, Container>* = nullptr>
// inline auto abs(const Container& x) {
//   return apply_vector_unary<Container>::apply(
//       x, [&](const auto& v) { return v.array().abs(); });
// }

namespace internal {
/**
 * Return the absolute value of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return absolute value of the argument
 */
template <typename V>
inline V complex_abs(const std::complex<V>& z) {
  return hypot(z.real(), z.imag());
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
