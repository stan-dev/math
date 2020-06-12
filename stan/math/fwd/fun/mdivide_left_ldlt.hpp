#ifndef STAN_MATH_FWD_FUN_MDIVIDE_LEFT_LDLT_HPP
#define STAN_MATH_FWD_FUN_MDIVIDE_LEFT_LDLT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/LDLT_factor.hpp>
#include <stan/math/prim/fun/mdivide_left_ldlt.hpp>
#include <stan/math/fwd/fun/to_fvar.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 *
 * @tparam R1 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C1 number of columns in the LDLT_factor, can be Eigen::Dynamic
 * @tparam EigMat type of the right-hand side matrix or vector
 *
 * @param A LDLT_factor
 * @param b right-hand side matrix or vector
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <int R1, int C1, typename EigMat,
          require_eigen_vt<is_fvar, EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, R1, EigMat::ColsAtCompileTime>
mdivide_left_ldlt(const LDLT_factor<double, R1, C1>& A, const EigMat& b) {
  using T = typename value_type_t<EigMat>::Scalar;
  constexpr int R2 = EigMat::RowsAtCompileTime;
  constexpr int C2 = EigMat::ColsAtCompileTime;
  check_multiplicable("mdivide_left_ldlt", "A", A, "b", b);

  const Eigen::Ref<const plain_type_t<EigMat>>& b_ref = b;
  Eigen::Matrix<T, R2, C2> b_val(b.rows(), b.cols());
  Eigen::Matrix<T, R2, C2> b_der(b.rows(), b.cols());
  for (int j = 0; j < b.cols(); j++) {
    for (int i = 0; i < b.rows(); i++) {
      b_val.coeffRef(i, j) = b_ref.coeff(i, j).val_;
      b_der.coeffRef(i, j) = b_ref.coeff(i, j).d_;
    }
  }

  return to_fvar(mdivide_left_ldlt(A, b_val), mdivide_left_ldlt(A, b_der));
}

}  // namespace math
}  // namespace stan
#endif
