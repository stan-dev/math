#ifndef STAN_MATH_PRIM_MAT_FUN_SUBTRACT_HPP
#define STAN_MATH_PRIM_MAT_FUN_SUBTRACT_HPP

#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/check_matching_dims.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return the result of subtracting the second specified matrix
 * from the first specified matrix.  The return scalar type is the
 * promotion of the input types.
 *
 * @tparam T1 Scalar type of first matrix.
 * @tparam T2 Scalar type of second matrix.
 * @tparam R Row type of matrices.
 * @tparam C Column type of matrices.
 * @param m1 First matrix.
 * @param m2 Second matrix.
 * @return Difference between first matrix and second matrix.
 */
template <typename Mat1, typename Mat2,
          typename = enable_if_all_eigen<Mat1, Mat2>>
inline auto subtract(const Mat1& m1, const Mat2& m2) {
  check_matching_dims("subtract", "m1", m1, "m2", m2);
  return m1 - m2;
}

template <typename Mat, typename Arith, typename = enable_if_eigen<Mat>,
          typename = enable_if_arithmetic<Arith>>
inline auto subtract(const Arith& c, const Mat& m) {
  return c - m.array();
}

template <typename Mat, typename Arith, typename = enable_if_arithmetic<Arith>,
          typename = enable_if_eigen<Mat>>
inline auto subtract(const Mat& m, const Arith& c) {
  return m.array() - c;
}

}  // namespace math
}  // namespace stan
#endif
