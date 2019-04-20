#ifndef STAN_MATH_PRIM_MAT_FUN_MDIVIDE_LEFT_TRI_HPP
#define STAN_MATH_PRIM_MAT_FUN_MDIVIDE_LEFT_TRI_HPP

#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/promote_common.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/lower_tri_inverse.hpp>
#include <stan/math/opencl/transpose.hpp>
#include <stan/math/opencl/copy.hpp>
#endif
namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b when A is triangular
 * @param A Triangular matrix.  Specify upper or lower with TriView
 * being Eigen::Upper or Eigen::Lower.
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <int TriView, typename T1, typename T2, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type,
                     R1, C2>
mdivide_left_tri(const Eigen::Matrix<T1, R1, C1> &A,
                 const Eigen::Matrix<T2, R2, C2> &b) {
  check_square("mdivide_left_tri", "A", A);
  check_multiplicable("mdivide_left_tri", "A", A, "b", b);
  return promote_common<Eigen::Matrix<T1, R1, C1>, Eigen::Matrix<T2, R1, C1> >(
             A)
      .template triangularView<TriView>()
      .solve(
          promote_common<Eigen::Matrix<T1, R2, C2>, Eigen::Matrix<T2, R2, C2> >(
              b));
}

/**
 * Returns the solution of the system Ax=b when A is triangular and b=I.
 * @param A Triangular matrix.  Specify upper or lower with TriView
 * being Eigen::Upper or Eigen::Lower.
 * @return x = A^-1 .
 * @throws std::domain_error if A is not square
 */
template <int TriView, typename T, int R1, int C1>
inline Eigen::Matrix<T, R1, C1> mdivide_left_tri(
    const Eigen::Matrix<T, R1, C1> &A) {
  check_square("mdivide_left_tri", "A", A);
  int n = A.rows();
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> b;
  b.setIdentity(n, n);
  A.template triangularView<TriView>().solveInPlace(b);
  return b;
}

#ifdef STAN_OPENCL
/**
 * Returns the solution of the system Ax=b when A is triangular
 * @param A Triangular matrix.  Specify upper or lower with TriView
 * being Eigen::Upper or Eigen::Lower.
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <int TriView, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<double, R1, C2>
mdivide_left_tri(const Eigen::Matrix<double, R1, C1> &A,
                 const Eigen::Matrix<double, R2, C2> &b) {
  check_square("mdivide_left_tri", "A", A);
  check_multiplicable("mdivide_left_tri", "A", A, "b", b);
  if (A.rows() >= 
      opencl_context.tuning_opts().lower_tri_inverse_size_worth_transfer) {
    matrix_cl A_cl(A);
    matrix_cl b_cl(b);
    matrix_cl A_inv_cl(A.rows(), A.cols());
    if (TriView == Eigen::Lower) {
      A_inv_cl = lower_triangular_inverse(A_cl);
    } else {
      A_inv_cl = transpose(lower_triangular_inverse(transpose(A_cl)));
    }    
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C(A.rows(),
                                                                 A.cols());
    auto C_cl = A_inv_cl * b_cl;
    copy(C, C_cl);  // NOLINT
    return C;
  } else {
    return promote_common<Eigen::Matrix<double, R1, C1>, Eigen::Matrix<double, R1, C1> >(
             A)
      .template triangularView<TriView>()
      .solve(
          promote_common<Eigen::Matrix<double, R2, C2>, Eigen::Matrix<double, R2, C2> >(
              b));  
  }  
}

/**
 * Returns the solution of the system Ax=b when A is triangular and b=I.
 * @param A Triangular matrix.  Specify upper or lower with TriView
 * being Eigen::Upper or Eigen::Lower.
 * @return x = A^-1 .
 * @note Because OpenCL only works on doubles there are two
 * <code>mdivide_left_tri(A)</code> functions. One that works on doubles
 * (this one) and another that works on all other types.
 * @throws std::domain_error if A is not square
 */
template <int TriView, int R1, int C1>
inline Eigen::Matrix<double, R1, C1> mdivide_left_tri(
    const Eigen::Matrix<double, R1, C1> &A) {
  check_square("mdivide_left_tri", "A", A);
  if (A.rows() >= 
      opencl_context.tuning_opts().lower_tri_inverse_size_worth_transfer) {
    matrix_cl A_cl(A);
    matrix_cl A_inv_cl(A.rows(), A.cols());
    if (TriView == Eigen::Lower) {
      A_inv_cl = lower_triangular_inverse(A_cl);
    } else {
      //TODO: Should we create an upper_triangular_inverse instea?
      // in terms of performance this is fine
      A_inv_cl = transpose(lower_triangular_inverse(transpose(A_cl)));
    }    
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A_inv(A.rows(),
                                                                 A.cols());
    copy(A_inv, A_inv_cl);  // NOLINT
    return A_inv;
  } else {
    int n = A.rows();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b;
    b.setIdentity(n, n);
    A.template triangularView<TriView>().solveInPlace(b);
    return b;
  }  
}
#endif

}  // namespace math
}  // namespace stan
#endif
