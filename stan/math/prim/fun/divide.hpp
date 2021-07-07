#ifndef STAN_MATH_PRIM_FUN_DIVIDE_HPP
#define STAN_MATH_PRIM_FUN_DIVIDE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cstddef>
#include <cstdlib>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Return the division of the first scalar by
 * the second scalar.
 * @param[in] x Specified scalar.
 * @param[in] y Specified scalar.
 * @return Scalar divided by the scalar.
 */
template <typename Scal1, typename Scal2,
          require_all_stan_scalar_t<Scal1, Scal2>* = nullptr>
inline return_type_t<Scal1, Scal2> divide(const Scal1& x, const Scal2& y) {
  return x / y;
}

inline int divide(int x, int y) {
  if (unlikely(y == 0)) {
    throw_domain_error("divide", "denominator is", y, "");
  }
  return x / y;
}

/**
 * Return matrix divided by scalar.
 *
 * @tparam Mat type of the matrix or expression
 * @tparam Scal type of the scalar
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename Mat, typename Scal, require_eigen_t<Mat>* = nullptr,
          require_stan_scalar_t<Scal>* = nullptr,
          require_all_not_var_t<scalar_type_t<Mat>, Scal>* = nullptr>
inline auto divide(const Mat& m, Scal c) {
  return m / c;
}

namespace internal {
template <typename T>
using is_fvar_or_arithmetic
    = bool_constant<std::is_arithmetic<scalar_type_t<T>>::value
                    || is_fvar<scalar_type_t<T>>::value>;
}
/**
 * Return matrix divided by matrix.
 *
 * @tparam Mat1 A type inheriting from `Eigen::EigenBase`
 * @tparam Mat2 A type inheriting from `Eigen::EigenBase`
 * @param[in] m specified matrix or expression
 * @param[in] c specified matrix or expression
 * @return matrix divided elementwise by `c`
 */
template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2>* = nullptr,
          require_all_t<internal::is_fvar_or_arithmetic<Mat1>,
                        internal::is_fvar_or_arithmetic<Mat2>>* = nullptr>
inline auto divide(const Mat1& m, const Mat2& c) {
  return (m.array() / c.array()).matrix();
}

/**
 * Return scalar divided by matrix.
 *
 * @tparam Mat A type inheriting from `Eigen::EigenBase`
 * @param[in] m specified matrix or expression
 * @param[in] c scalar double
 * @return matrix divided elementwise by `c`
 */
template <typename Mat, require_eigen_vt<std::is_arithmetic, Mat>* = nullptr>
inline auto divide(double c, const Mat& m) {
  return (c / m.array()).matrix();
}

/**
 * Return scalar divided by matrix.
 *
 * @tparam Scalar A scalar of Arithmetic of `fvar` type.
 * @tparam Mat A type inheriting from `Eigen::EigenBase`
 * @param[in] c scalar
 * @param[in] m specified matrix or expression
 * @return matrix divided elementwise by `c`
 */
template <typename Scalar, typename Mat, require_eigen_t<Mat>* = nullptr,
          require_t<bool_constant<std::is_arithmetic<Scalar>::value
                                  || is_fvar<Scalar>::value>>* = nullptr,
          require_any_st_fvar<Scalar, Mat>* = nullptr>
inline auto divide(Scalar c, const Mat& m) {
  return (c / m.array()).matrix();
}

}  // namespace math
}  // namespace stan

#endif
