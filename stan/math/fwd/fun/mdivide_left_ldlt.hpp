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
 * @tparam T type of matrix for the LDLT_factor
 * @tparam EigMat type of the right-hand side matrix or vector
 * @param A LDLT_factor
 * @param b right-hand side matrix or vector
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <typename T, typename EigMat,
          require_eigen_vt<std::is_arithmetic, T>* = nullptr,
          require_eigen_vt<is_fvar, EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic,
                     EigMat::ColsAtCompileTime>
mdivide_left_ldlt(LDLT_factor<T>& A, const EigMat& b) {
  using EigMatValueScalar = typename value_type_t<EigMat>::Scalar;
  constexpr int R2 = EigMat::RowsAtCompileTime;
  constexpr int C2 = EigMat::ColsAtCompileTime;
  check_multiplicable("mdivide_left_ldlt", "A", A.matrix(), "b", b);

  const auto& b_ref = to_ref(b);
  Eigen::Matrix<EigMatValueScalar, R2, C2> b_val(b.rows(), b.cols());
  Eigen::Matrix<EigMatValueScalar, R2, C2> b_der(b.rows(), b.cols());
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
