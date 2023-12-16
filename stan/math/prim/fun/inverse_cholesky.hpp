#ifndef STAN_MATH_PRIM_FUN_INVERSE_CHOLESKY_HPP
#define STAN_MATH_PRIM_FUN_INVERSE_CHOLESKY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/mdivide_left_tri.hpp>
#include <stan/math/prim/fun/inv_square.hpp>

namespace stan {
namespace math {

/**
 * Returns the inverse of the matrix whose Cholesky factor is L
 *
 * @tparam T type of elements in the matrix
 * @param L Matrix that is a Cholesky factor.
 * @return The matrix inverse of L * L'
 * @throw std::domain_error If the input matrix is not square or
 *  lower triangular
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_eigen_vt<is_var, EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::ColsAtCompileTime>
inverse_cholesky(const EigMat& L) {
  const Eigen::Ref<const plain_type_t<EigMat>>& L_ref = L;
  check_square("inverse_cholesky", "L", L_ref);
  check_lower_triangular("inverse_cholesky", "L", L_ref);
  int K = L.rows();
  using T_result = plain_type_t<EigMat>;
  
  if (K == 0) {
    return {};
  }
  if (K == 1) {
    T_result X(1, 1);
    X.coeffRef(0) = 1 / L_ref.coeff(0, 0);
    return X;
  }
  return mdivide_left_tri<Eigen::Lower>(L_ref);
}

}  // namespace math
}  // namespace stan

#endif
