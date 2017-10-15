#ifndef STAN_MATH_PRIM_MAT_FUN_CHOL2INV_HPP
#define STAN_MATH_PRIM_MAT_FUN_CHOL2INV_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/scal/fun/inv_square.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_lower_triangular.hpp>

namespace stan {
  namespace math {

    /**
     * Returns the inverse of the matrix whose Cholesky factor is L
     * @param L Matrix that is a Cholesky factor.
     * @return The matrix inverse of L * L'
     * @throw std::domain_error If the input matrix is not square or
     *  lower triangular
     */
    inline matrix_d
    chol2inv(const matrix_d& L) {
      check_square("chol2inv", "L", L);
      check_lower_triangular("chol2inv", "L", L);
      int K = L.rows();
      if (K == 0)
        return L;
      if (K == 1) {
        matrix_d X(1, 1);
        X.coeffRef(0) = inv_square(L.coeff(0));
        return X;
      }
      matrix_d L_inv = L .template triangularView<Eigen::Lower>()
                         .solve(matrix_d::Identity(K,K));
      matrix_d X(K, K);
      for (int k = 0; k < K; ++k) {
        X.coeffRef(k, k) = L_inv.col(k).tail(K - k).squaredNorm();
        for (int j = k + 1; j < K; ++j) {
          int Kmj = K - j;
          X.coeffRef(k, j) = X.coeffRef(j, k) =
            L_inv.col(k).tail(Kmj).dot(L_inv.col(j).tail(Kmj));
        }
      }
      return X;
    }

  }
}
#endif
