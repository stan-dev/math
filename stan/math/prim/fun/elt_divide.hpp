#ifndef STAN_MATH_PRIM_FUN_ELT_DIVIDE_HPP
#define STAN_MATH_PRIM_FUN_ELT_DIVIDE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/divide.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise division of the specified matrices.
 *
 * @tparam Mat1 type of the first matrix or expression
 * @tparam Mat2 type of the second matrix or expression
 *
 * @param m1 First matrix or expression
 * @param m2 Second matrix or expression
 * @return Elementwise division of matrices.
 */
template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2>* = nullptr,
          require_all_not_st_var<Mat1, Mat2>* = nullptr>
auto elt_divide(const Mat1& m1, const Mat2& m2) {
  check_matching_dims("elt_divide", "m1", m1, "m2", m2);
#ifdef USE_STANC3
  return (m1.array() / m2.array()).matrix();
#else
  return (m1.array() / m2.array()).matrix().eval();
#endif
}

/**
 * Return the elementwise division of the specified matrix
 * by the specified scalar.
 *
 * @tparam Mat type of the matrix or expression
 * @tparam Scal type of the scalar
 *
 * @param m matrix or expression
 * @param s scalar
 * @return Elementwise division of a scalar by matrix.
 */
#ifdef USE_STANC3
template <typename Mat, typename Scal, require_matrix_t<Mat>* = nullptr,
          require_stan_scalar_t<Scal>* = nullptr>
auto elt_divide(const Mat& m, Scal s) {
  return divide(m, s);
}
#else
template <typename Scal, typename Mat, typename = require_stan_scalar_t<Scal>,
          typename = require_eigen_t<Mat> >
auto elt_divide(const Mat& m, Scal s) {
  return (m.array() / s).matrix().eval();
}
#endif

/**
 * Return the elementwise division of the specified scalar
 * by the specified matrix.
 *
 * @tparam Scal type of the scalar
 * @tparam Mat type of the matrix or expression
 *
 * @param s scalar
 * @param m matrix or expression
 * @return Elementwise division of a scalar by matrix.
 */
template <typename Scal, typename Mat, require_stan_scalar_t<Scal>* = nullptr,
          require_eigen_t<Mat>* = nullptr>
auto elt_divide(Scal s, const Mat& m) {
#ifdef USE_STANC3
  return (s / m.array()).matrix();
#else
  return (s / m.array()).matrix().eval();
#endif
}

template <typename Scal1, typename Scal2,
          require_all_stan_scalar_t<Scal1, Scal2>* = nullptr>
auto elt_divide(Scal1 s1, Scal2 s2) {
  return s1 / s2;
}

#ifndef USE_STANC3
/**
 * Return the elementwise division of the specified matrices.
 *
 * @tparam T1 Type of scalars in first matrix.
 * @tparam T2 Type of scalars in second matrix.
 * @tparam R Row type of both matrices.
 * @tparam C Column type of both matrices.
 * @param m1 First matrix
 * @param m2 Second matrix
 * @return Elementwise division of matrices.
 */
template <typename T1, typename T2, int R, int C>
Eigen::Matrix<return_type_t<T1, T2>, R, C> elt_divide(
    const Eigen::Matrix<T1, R, C>& m1, const Eigen::Matrix<T2, R, C>& m2) {
  check_matching_dims("elt_divide", "m1", m1, "m2", m2);

  return m1.array() / m2.array();
}

/**
 * Return the elementwise division of the specified matrix
 * by the specified scalar.
 *
 * @tparam T1 Type of scalars in the matrix.
 * @tparam T2 Type of the scalar.
 * @tparam R Row type of the matrix.
 * @tparam C Column type of the matrix.
 * @param m matrix
 * @param s scalar
 * @return Elementwise division of a scalar by matrix.
 */
template <typename T1, typename T2, int R, int C,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
Eigen::Matrix<return_type_t<T1, T2>, R, C> elt_divide(
    const Eigen::Matrix<T1, R, C>& m, T2 s) {
  return m.array() / s;
}

/**
 * Return the elementwise division of the specified scalar
 * by the specified matrix.
 *
 * @tparam T1 Type of the scalar.
 * @tparam T2 Type of scalars in the matrix.
 * @tparam R Row type of the matrix.
 * @tparam C Column type of the matrix.
 * @param s scalar
 * @param m matrix
 * @return Elementwise division of a scalar by matrix.
 */
template <typename T1, typename T2, int R, int C,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
Eigen::Matrix<return_type_t<T1, T2>, R, C> elt_divide(
    T1 s, const Eigen::Matrix<T2, R, C>& m) {
  return s / m.array();
}
#endif

}  // namespace math
}  // namespace stan

#endif
