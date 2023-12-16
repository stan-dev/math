#ifndef STAN_MATH_REV_FUN_INVERSE_CHOLESKY_HPP
#define STAN_MATH_REV_FUN_INVERSE_CHOLESKY_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inverse_cholesky.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <memory>

namespace stan {
namespace math {

// namespace internal {
// template <typename LMat, typename LAMat>
// inline void initialize_return(LMat& L, const LAMat& L_A, vari*& dummy) {
//   for (Eigen::Index j = 0; j < L_A.rows(); ++j) {
//     for (Eigen::Index i = 0; i < L_A.rows(); ++i) {
//       if (j > i) {
//         L.coeffRef(i, j) = dummy;
//       } else {
//         L.coeffRef(i, j) = new vari(L_A.coeffRef(i, j), false);
//       }
//     }
//   }
// }

// template <typename T1, typename T2>
// inline auto cholesky_inverse_lambda(T1& arena_A, T2& inv_A) {
//   return [arena_A, inv_A]() mutable {
//     arena_A.adj()
//         += 
//   };
// }

// }


/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 *
 * @tparam T type of B
 * @param A LDLT_factor
 * @param B Right hand side matrix or vector.
 * @return x = A^-1 B, solution of the linear system.
 * @throws std::domain_error if rows of B don't match the size of A.
 */
template <typename T, require_eigen_vt<is_var, T>*  = nullptr>
inline auto inverse_cholesky(const T& A) {
   using ret_type = return_var_matrix_t<T>;

   if (unlikely(A.size() == 0)) {
    return A;
  }

    arena_t<promote_scalar_t<var, T>> arena_A = A;
    auto arena_A_val = to_arena(arena_A.val());
    arena_t<ret_type> res = inverse_cholesky(arena_A_val);

    reverse_pass_callback([arena_A, arena_A_val, res]() mutable {
  
      // arena_A.adj() -= (res_val * res.val().transpose().eval())
      //                      .template triangularView<Eigen::Lower>();
     promote_scalar_t<double, T> adjB
          = arena_A_val.template triangularView<Eigen::Lower>().transpose().solve(
              res.adj());
      arena_A.adj() -= (adjB * res.val().transpose().eval())
                           .template triangularView<Eigen::Lower>();
    });

    return plain_type_t<T>(res);
}

// template <Eigen::UpLoType TriView, typename T, require_eigen_vt<is_var, T>* = nullptr>
// inline auto inverse_cholesky(const T &A) {
//   using ret_val_type
//       = Eigen::Matrix<double, T::RowsAtCompileTime, T::ColsAtCompileTime>;
//   using ret_type = promote_var_matrix_t<ret_val_type, T, T>;

//   if (A.matrix().size() == 0) {
//     return ret_type(ret_val_type(0, 0));
//   }

//     arena_t<T> arena_A = A;
//     arena_t<promote_scalar_t<double, T>> res_val = arena_A.val().inverse();
//     arena_t<ret_type> res = res_val;

//     reverse_pass_callback([arena_A, res_val, res]() mutable {
//      arena_A.adj() -= res_val.transpose() * res.adj_op() * res_val.transpose();
//     });

//     return ret_type(res);
// }

}  // namespace math
}  // namespace stan
#endif
