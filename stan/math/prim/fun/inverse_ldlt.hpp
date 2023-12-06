#ifndef STAN_MATH_PRIM_FUN_INVERSE_LDLT_HPP
#define STAN_MATH_PRIM_FUN_INVERSE_LDLT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/identity_matrix.hpp>
#include <stan/math/prim/fun/LDLT_factor.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A and b=I.
 *
 * @tparam T type of the matrix
 *
 * @param A LDLT_factor
 * @return x = A^-1 .
 * @throws std::domain_error if A is not square
 */
template <typename T, require_eigen_t<T> * = nullptr>
inline plain_type_t<T> inverse_ldlt(LDLT_factor<T>& A) {

  if (A.matrix().cols() == 0) {
    return {0, 0};
  }
  const Eigen::Ref<const plain_type_t<T>>& A_ref = A.matrix(); 
  const int n = A_ref.rows();
  
  plain_type_t<T> b = plain_type_t<T>::Identity(n, n);
  
  A_ref.ldlt().solveInPlace(b);
  
  return b;
}

}  // namespace math
}  // namespace stan

#endif
