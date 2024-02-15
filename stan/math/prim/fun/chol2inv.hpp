#ifndef STAN_MATH_PRIM_FUN_CHOL2INV_HPP
#define STAN_MATH_PRIM_FUN_CHOL2INV_HPP

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
template <typename T, require_eigen_t<T>* = nullptr>
plain_type_t<T> chol2inv(const T& L) {
  const Eigen::Ref<const plain_type_t<T>>& L_ref = L;
  check_square("chol2inv", "L", L_ref);
  check_lower_triangular("chol2inv", "L", L_ref);
  int K = L.rows();
  using T_result = plain_type_t<T>;
  if (K == 0) {
    return L_ref;
  }
  if (K == 1) {
    T_result X(1, 1);
    X.coeffRef(0) = inv_square(L_ref.coeff(0, 0));
    return X;
  }
  T_result L_inv = mdivide_left_tri<Eigen::Lower>(L_ref);
  T_result X(K, K);
  for (int k = 0; k < K; ++k) {
    X.coeffRef(k, k) = dot_self(L_inv.col(k).tail(K - k).eval());
    for (int j = k + 1; j < K; ++j) {
      int Kmj = K - j;
      X.coeffRef(k, j) = X.coeffRef(j, k) = dot_product(
          L_inv.col(k).tail(Kmj).eval(), L_inv.col(j).tail(Kmj).eval());
    }
  }
  return X;
}

}  // namespace math
}  // namespace stan

#endif
