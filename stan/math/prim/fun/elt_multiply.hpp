#ifndef STAN_MATH_PRIM_FUN_ELT_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_ELT_MULTIPLY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise multiplication of the specified
 * matrices.
 *
 * @tparam Mat1 type of the first matrix or expression
 * @tparam Mat2 type of the second matrix or expression
 *
 * @param m1 First matrix or expression
 * @param m2 Second matrix or expression
 * @return Elementwise product of matrices.
 */
template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2>* = nullptr,
          require_all_not_st_var<Mat1, Mat2>* = nullptr>
auto elt_multiply(const Mat1& m1, const Mat2& m2) {
  check_matching_dims("elt_multiply", "m1", m1, "m2", m2);
  return m1.cwiseProduct(m2);
}

/**
 * Return the elementwise multiplication of the specified
 * matrix and a scalar.
 *
 * @tparam Mat type of the matrix or expression
 * @tparam Scal type of the scalar
 *
 * @param m matrix or expression
 * @param s scalar
 * @return Elementwise product of a matrix and a scalar.
 */
template <typename Mat, typename Scal,
          require_eigen_t<Mat>* = nullptr,
          require_stan_scalar_t<Scal>* = nullptr>
auto elt_multiply(const Mat& m, const Scal& s) {
  return m * s;
}

/**
 * Return the elementwise multiplication of the specified
 * matrix and a scalar.
 *
 * @tparam Mat1 type of the matrix or expression
 * @tparam Mat2 type of the scalar
 *
 * @param m1 matrix or expression
 * @param m2 scalar
 * @return Elementwise product of a matrix and a scalar.
 */
template <typename Mat, typename Scal,
          require_eigen_t<Mat>* = nullptr,
          require_stan_scalar_t<Scal>* = nullptr>
auto elt_multiply(const Scal& s, const Mat& m) {
  return s * m;
}

/**
 * Return the elementwise multiplication of the specified
 * scalars.
 *
 * @tparam Mat1 type of the first scalar
 * @tparam Mat2 type of the second scalar
 *
 * @param a First scalar
 * @param b Second scalar
 * @return Product of scalars
 */
template <typename Scalar1, typename Scalar2,
          require_all_stan_scalar_t<Scalar1, Scalar2>* = nullptr>
auto elt_multiply(const Scalar1& a, const Scalar2& b) {
  return a * b;
}

}  // namespace math
}  // namespace stan

#endif
