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
  return (m1.array() / m2.array()).matrix();
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
template <typename Mat, typename Scal, require_matrix_t<Mat>* = nullptr,
          require_stan_scalar_t<Scal>* = nullptr>
auto elt_divide(const Mat& m, Scal s) {
  return divide(m, s);
}

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
  return (s / m.array()).matrix();
}

template <typename Scal1, typename Scal2,
          require_all_stan_scalar_t<Scal1, Scal2>* = nullptr>
auto elt_divide(Scal1 s1, Scal2 s2) {
  return s1 / s2;
}

}  // namespace math
}  // namespace stan

#endif
