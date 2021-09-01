#ifndef STAN_MATH_PRIM_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_PRIM_FUN_VALUE_OF_REC_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <complex>
#include <cstddef>
#include <vector>

namespace stan {
namespace math {


/**
 * Return the specified argument.
 *
 * <p>See <code>value_of(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @param x Specified value.
 * @return Specified value.
 */
inline double value_of_rec(double x) { return x; }

/**
 * Recursively apply value-of to the parts of the argument.
 *
 * @tparam T value type of argument
 * @param[in] x argument
 * @return real complex value of argument
 */
template <typename T>
inline std::complex<double> value_of_rec(const std::complex<T>& x) {
  return {value_of_rec(x.real()), value_of_rec(x.imag())};
}

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of_rec(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @param x Specified std::vector.
 * @return Specified std::vector.
 */
template <typename T, require_std_vector_t<T>* = nullptr,
          require_st_same<double, T>* = nullptr>
inline T value_of_rec(T&& x) {
  return std::forward<T>(x);
}

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of_rec(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @tparam T Type of matrix.
 * @param x Specified matrix.
 * @return Specified matrix.
 */
template <typename T, require_eigen_t<T>* = nullptr,
  require_st_same<T, double>* = nullptr>
inline T value_of_rec(T&& x) {
  return std::forward<T>(x);
}

/**
 * Convert a matrix of type T to a matrix of doubles.
 *
 * T must implement value_of_rec. See
 * test/unit/math/fwd/fun/value_of_test.cpp for fvar and var usage.
 *
 * @tparam T Type of matrix
 * @param[in] M Matrix to be converted
 * @return Matrix of values
 **/
template <typename T, require_not_st_same<T, double>* = nullptr,
          require_eigen_t<T>* = nullptr>
inline auto value_of_rec(T&& M) {
  return make_holder(
      [](auto& m) {
        return m.unaryExpr([](auto x) { return value_of_rec(x); });
      },
      std::forward<T>(M));
}

/**
 * Convert a std::vector of type T to an std::vector of T with a double scalar
 *type.
 *
 * T must implement value_of_rec. See
 * test/math/fwd/fun/value_of_rec.cpp for fvar and var usage.
 *
 * @tparam T type in std::vector
 * @param[in] x std::vector to be extract values from.
 * @return std::vector of values
 **/
template <typename T, require_std_vector_t<T>* = nullptr, require_not_st_same<double, T>* = nullptr>
inline auto value_of_rec(const T& x) {
  promote_scalar_t<double, T> result(x.size());
  std::transform(x.begin(), x.end(), result.begin(),
                 [](auto&& xx) { return value_of_rec(xx); });
  return result;
}

}  // namespace math
}  // namespace stan

#endif
