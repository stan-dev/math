#ifndef STAN_MATH_FWD_FUN_POW_HPP
#define STAN_MATH_FWD_FUN_POW_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/sqrt.hpp>
#include <stan/math/fwd/fun/inv.hpp>
#include <stan/math/fwd/fun/inv_sqrt.hpp>
#include <stan/math/fwd/fun/inv_square.hpp>
#include <stan/math/prim/fun/pow.hpp>
#include <cmath>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> pow(const fvar<T>& x1, const fvar<T>& x2) {
  using std::log;
  using std::pow;
  T pow_x1_x2(pow(x1.val_, x2.val_));
  return fvar<T>(pow_x1_x2, (x2.d_ * log(x1.val_) + x2.val_ * x1.d_ / x1.val_)
                                * pow_x1_x2);
}

template <typename T, typename U, typename = require_arithmetic_t<U>>
inline fvar<T> pow(U x1, const fvar<T>& x2) {
  using std::log;
  using std::pow;
  T u = pow(x1, x2.val_);
  return fvar<T>(u, x2.d_ * log(x1) * u);
}

template <typename T, typename U, typename = require_arithmetic_t<U>>
inline fvar<T> pow(const fvar<T>& x1, U x2) {
  using std::pow;
  using std::sqrt;
  if (x2 == -2) {
    return inv_square(x1);
  }
  if (x2 == -1) {
    return inv(x1);
  }
  if (x2 == -0.5) {
    return inv_sqrt(x1);
  }
  if (x2 == 0.5) {
    return sqrt(x1);
  }
  if (x2 == 1.0) {
    return x1;
  }
  if (x2 == 2.0) {
    return square(x1);
  }
  return fvar<T>(pow(x1.val_, x2), x1.d_ * x2 * pow(x1.val_, x2 - 1));
}

// must uniquely match all pairs of:
//    { complex<fvar<V>>, complex<T>, fvar<V>, T }
// with at least one fvar<V> and at least one complex, where T is arithmetic:
// 1) complex<fvar<V>>, complex<fvar<V>>
// 2) complex<fvar<V>>, complex<T>
// 3) complex<fvar<V>>, fvar<V>
// 4) complex<fvar<V>>, T
// 5) complex<T>, complex<fvar<V>>
// 6) complex<T>, fvar<V>
// 7) fvar<V>, complex<fvar<V>>
// 8) fvar<V>, complex<T>
// 9) T, complex<fvar<V>>

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename V>
inline std::complex<fvar<V>> pow(const std::complex<fvar<V>>& x,
                                 const std::complex<fvar<V>>& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam V autodiff value type
 * @tparam T arithmetic type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename V, typename T, typename = require_arithmetic_t<T>>
inline std::complex<fvar<V>> pow(const std::complex<fvar<V>>& x,
                                 const std::complex<T>& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam V autodiff value type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename V>
inline std::complex<fvar<V>> pow(const std::complex<fvar<V>>& x,
                                 const fvar<V>& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam V autodiff value type
 * @tparam T arithmetic type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename V, typename T, typename = require_arithmetic_t<T>>
inline std::complex<fvar<V>> pow(const std::complex<fvar<V>>& x, const T& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam V autodiff value type
 * @tparam T arithmetic type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename V, typename T, typename = require_arithmetic_t<T>>
inline std::complex<fvar<V>> pow(const std::complex<T>& x,
                                 const std::complex<fvar<V>>& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam V autodiff value type
 * @tparam T arithmetic type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename V, typename T, typename = require_arithmetic_t<T>>
inline std::complex<fvar<V>> pow(const std::complex<T>& x, const fvar<V>& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam V autodiff value type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename V>
inline std::complex<fvar<V>> pow(const fvar<V>& x,
                                 const std::complex<fvar<V>>& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam V autodiff value type
 * @tparam T arithmetic type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename V, typename T, typename = require_arithmetic_t<T>>
inline std::complex<fvar<V>> pow(const fvar<V>& x, const std::complex<T>& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam V autodiff value type
 * @tparam T real type (`fvar<V>` or arithmetic)
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename T, typename V, typename = require_arithmetic_t<T>>
inline std::complex<fvar<V>> pow(T x, const std::complex<fvar<V>>& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * Note: this overload is required because gcc still provides the
 * C++99 template function `pow(complex<T>, int)`, which introduces
 * an ambiguity.
 *
 * @tparam T autodiff value type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename T>
inline std::complex<fvar<T>> pow(const std::complex<fvar<T>>& x, int y) {
  return internal::complex_pow(x, y);
}

}  // namespace math
}  // namespace stan
#endif
