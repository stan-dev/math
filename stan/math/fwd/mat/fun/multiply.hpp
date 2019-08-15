#ifndef STAN_MATH_FWD_MAT_FUN_MULTIPLY_HPP
#define STAN_MATH_FWD_MAT_FUN_MULTIPLY_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/mat/fun/typedefs.hpp>
#include <stan/math/fwd/mat/fun/dot_product.hpp>

namespace stan {
namespace math {

  /**
   * Return specified matrix multiplied by specified scalar.
   * @tparam R Row type for matrix.
   * @tparam C Column type for matrix.
   * @param m Matrix.
   * @param c Scalar.
   * @return Product of matrix and scalar.
   */
  template <typename T1, typename T2, enable_if_any_fvar<T1, T2>* = nullptr>
  inline auto multiply(T1 m, T2 c) {
    return c * m;
  }

  /**
   * Return the product of the specified matrices.  The number of
   * columns in the first matrix must be the same as the number of rows
   * in the second matrix.
   * @param m1 First matrix.
   * @param m2 Second matrix.
   * @return The product of the first and second matrices.
   * @throw std::domain_error if the number of columns of m1 does not match
   *   the number of rows of m2.
   */
  template <typename T1, typename T2, enable_if_all_eigen<T1, T2>* = nullptr,
            enable_if_all_fvar<scalar_type_t<T1>, scalar_type_t<T2>>* = nullptr,
            enable_if_not_dot_product<T1, T2>* = nullptr>
  inline auto multiply(const T1& m1, const T2& m2) {
    check_multiplicable("multiply", "m1", m1, "m2", m2);
    return m1 * m2;
  }

  /**
   * Return the scalar product of the specified row vector and
   * specified column vector.  The return is the same as the dot
   * product.  The two vectors must be the same size.
   * @param rv Row vector.
   * @param v Column vector.
   * @return Scalar result of multiplying row vector by column vector.
   * @throw std::domain_error if rv and v are not the same size.
   */
  template <typename T1, typename T2, enable_if_all_fvar<scalar_type_t<T1>, scalar_type_t<T2>>* = nullptr,
            enable_if_dot_product<T1, T2>* = nullptr>
  inline auto multiply(const T1& rv, const T2& v) {
    check_matching_sizes("multiply", "rv", rv, "v", v);
    return rv.dot(v);
  }



}  // namespace math
}  // namespace stan
#endif
