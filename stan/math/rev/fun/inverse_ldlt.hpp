#ifndef STAN_MATH_REV_FUN_INVERSE_LDLT_HPP
#define STAN_MATH_REV_FUN_INVERSE_LDLT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/LDLT_factor.hpp>
#include <stan/math/prim/fun/inverse_ldlt.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <memory>

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 *
 * @tparam T type of B
 * @param A LDLT_factor
 * @param B Right hand side matrix or vector.
 * @return x = A^-1 B, solution of the linear system.
 * @throws std::domain_error if rows of B don't match the size of A.
 */
template <typename T, require_var_matrix_t<T> * = nullptr>
inline auto inverse_ldlt(LDLT_factor<T>& A) {
  using ret_val_type
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  using ret_type = return_var_matrix_t<T>;

  if (A.matrix().size() == 0) {
    return ret_type(ret_val_type(0, 0));
  }
   const int n = A.matrix().rows();

    arena_t<promote_scalar_t<var, T>> arena_A = A.matrix();
    arena_t<ret_type> res = inverse_ldlt(A);
    // arena_t<ret_type> res = A.ldlt().solve(Eigen::MatrixXd::Identity(n, n));
    const auto* ldlt_ptr = make_chainable_ptr(A.ldlt());

    reverse_pass_callback([arena_A, ldlt_ptr, res]() mutable {
      arena_A.adj() -= ldlt_ptr->solve(res.adj()) * res.val_op().transpose();
    });

    return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif