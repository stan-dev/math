#ifndef STAN_MATH_PRIM_FUN_KRONECKER_HPP
#define STAN_MATH_PRIM_FUN_KRONECKER_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Return the Kronecker product of two matrices, the block matrix whose
 * block \f$(i, j)\f$ is \f$A_{ij} B\f$:
 * \f[
 * (A \otimes B)_{p i + k,\, q j + l} = A_{ij} B_{kl},
 * \f]
 * where \f$B\f$ is \f$p \times q\f$.
 *
 * Column \f$q j + l\f$ of the result, read as a \f$p \times m\f$ matrix, is
 * the outer product of column \f$l\f$ of \f$B\f$ and column \f$j\f$ of
 * \f$A\f$; the result is filled one such column at a time.
 *
 * @tparam T1 type of the first matrix
 * @tparam T2 type of the second matrix
 * @param A first matrix (m x n)
 * @param B second matrix (p x q)
 * @return Kronecker product of `A` and `B` (mp x nq)
 */
template <typename T1, typename T2,
          require_all_eigen_matrix_base_vt<is_stan_scalar, T1, T2>* = nullptr,
          require_not_return_type_t<is_var, T1, T2>* = nullptr>
inline Eigen::Matrix<return_type_t<T1, T2>, Eigen::Dynamic, Eigen::Dynamic>
kronecker(T1&& A, T2&& B) {
  decltype(auto) A_ref = to_ref(std::forward<T1>(A));
  decltype(auto) B_ref = to_ref(std::forward<T2>(B));
  const Eigen::Index m = A_ref.rows();
  const Eigen::Index p = B_ref.rows();
  const Eigen::Index q = B_ref.cols();
  Eigen::Matrix<return_type_t<T1, T2>, Eigen::Dynamic, Eigen::Dynamic> res(
      m * p, A_ref.cols() * q);
  for (Eigen::Index c = 0; c < res.cols(); ++c) {
    res.col(c).reshaped(p, m).noalias()
        = B_ref.col(c % q) * A_ref.col(c / q).transpose();
  }
  return res;
}

}  // namespace math
}  // namespace stan

#endif
