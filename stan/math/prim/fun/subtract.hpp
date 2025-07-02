#ifndef STAN_MATH_PRIM_FUN_SUBTRACT_HPP
#define STAN_MATH_PRIM_FUN_SUBTRACT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the result of subtracting the second scalar from the first
 * scalar.
 *
 * @tparam ScalarA type of the first scalar
 * @tparam ScalarB type of the second scalar
 * @param a first scalar
 * @param b second scalar
 * @return difference between first scalar and second scalar
 */
template <typename ScalarA, typename ScalarB,
          require_all_stan_scalar_t<ScalarA, ScalarB>* = nullptr,
          require_all_not_var_t<ScalarA, ScalarB>* = nullptr>
inline return_type_t<ScalarA, ScalarB> subtract(const ScalarA& a,
                                                const ScalarB& b) {
  return a - b;
}

/**
 * Return the result of subtracting the second specified matrix
 * from the first specified matrix.
 *
 * @tparam Mat1 type of the first matrix or expression
 * @tparam Mat2 type of the second matrix or expression
 *
 * @param m1 First matrix or expression.
 * @param m2 Second matrix or expression.
 * @return Difference between first matrix and second matrix.
 */
template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2>* = nullptr,
          require_all_not_st_var<Mat1, Mat2>* = nullptr>
inline auto subtract(Mat1&& m1, Mat2&& m2) {
  check_matching_dims("subtract", "m1", m1, "m2", m2);
  return make_holder([](auto&& m1_, auto&& m2_) {
    return std::forward<decltype(m1_)>(m1_)
           - std::forward<decltype(m2_)>(m2_);
  }, std::forward<Mat1>(m1), std::forward<Mat2>(m2));
}

/**
 * Return the result of subtracting the specified matrix from the specified
 * scalar.
 *
 * @tparam Scal type of the scalar
 * @tparam Mat type of the matrix or expression
 * @param c Scalar.
 * @param m Matrix or expression.
 * @return The scalar minus the matrix.
 */
template <typename Scal, typename Mat, require_stan_scalar_t<Scal>* = nullptr,
          require_eigen_t<Mat>* = nullptr,
          require_all_not_st_var<Mat, Scal>* = nullptr>
inline auto subtract(const Scal c, Mat&& m) {
  return make_holder([c](auto&& m_) {
    return (c - std::forward<decltype(m_)>(m_).array()).matrix();
  }, std::forward<Mat>(m));
}

/**
 * Return the result of subtracting the specified scalar from the specified
 * matrix.
 *
 * @tparam Mat type of the matrix or expression
 * @tparam Scal type of the scalar
 * @param m Matrix or expression.
 * @param c Scalar.
 * @return The matrix minus the scalar.
 */
template <typename Mat, typename Scal, require_eigen_t<Mat>* = nullptr,
          require_stan_scalar_t<Scal>* = nullptr,
          require_all_not_st_var<Scal, Mat>* = nullptr>
inline auto subtract(Mat&& m, const Scal c) {
  return make_holder([c](auto&& m_) {
    return (std::forward<decltype(m_)>(m_).array() - c).matrix();
  }, std::forward<Mat>(m));
}

}  // namespace math
}  // namespace stan

#endif
