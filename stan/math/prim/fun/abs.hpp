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
  
/**
 * Metaprogram calculating return type for applying `abs` to an
 * argument of the specified template type.  This struct defines the
 * default case, which will be used for arithmetic types, with
 * specializations dealing with other cases.
 *
 * @tparam T argument type
 */
template <typename T>
struct abs_return {  
  /**
   * Return type of `abs(T)`.
   */
  using type = T;
};

/**
 * Helper typedef for abs_return.  With this definition,
 * `typename abs_return<T>::type` can be abbreviated to `abs_return_t<T>`.
 *
 * @tparam T type of argument
 */ 
template <typename T>
using abs_return_t = typename abs_return<T>::type;

template <typename T>
struct abs_return<std::complex<T>> {
  using type = T;
};

template <typename T, int R, int C>
struct abs_return<Eigen::Matrix<T, R, C>> {
  using type = Eigen::Matrix<abs_return_t<T>, R, C>;
};

template <typename T>
struct abs_return<std::vector<T>> {
  using type = std::vector<abs_return_t<T>>;
};

/**
 * Return the absolute value of the specified arithmetic argument.
 *
 * @tparam T type of argument (must be arithmetic)
 * @param x argument
 * @return absolute value of argument
 */
template <typename T>
inline auto abs(T x) {
  return std::abs(x);
}

/**
 * Return the elementwise absolute value of the specified matrix or vector argument.
 *
 * @tparam T type of matrix elements
 * @param x argument
 * @return elementwise absolute value of argument
 */
template <typename T, int R, int C>
inline auto abs(const Eigen::Matrix<T, R, C>& x) {
  Eigen::Matrix<abs_return_t<T>, R, C> y(x.rows(), x.cols());
  for (int i = 0; i < x.size(); ++i)
    y(i) = abs(x(i));
  return y;
}

/**
 * Return the elementwise absolute value of the specified standard vector argument.
 *
 * @tparam T type of vector elements
 * @param x argument
 * @return elementwise absolute value of argument
 */
template <typename T>
inline auto abs(const std::vector<T>& x) {
  std::vector<abs_return_t<T>> y;
  y.reserve(x.size());
  for (const auto& xi : x)
    y.push_back(abs(xi));
  return y;
}
  
}  // namespace math
}  // namespace stan

#endif
