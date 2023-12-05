#ifndef STAN_MATH_REV_FUN_INVERSE_LDLT_HPP
#define STAN_MATH_REV_FUN_INVERSE_LDLT_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/inverse_ldlt.hpp>
#include <stan/math/rev/fun/LDLT_factor.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Reverse mode specialization of calculating the inverse of the matrix.
 *
 * @param A specified matrix
 * @return Inverse of the matrix (an empty matrix if the specified matrix has
 * size zero).
 * @throw std::invalid_argument if the matrix is not square.
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto inverse_ldlt(LDLT_factor<T>& A) {
  using ret_type = return_var_matrix_t<T>;
  if (unlikely(A.matrix().size() == 0)) {
    return ret_type(A.matrix());
  }

  arena_t<T> arena_A = A.matrix();
  arena_t<promote_scalar_t<double, T>> res_val = inverse_ldlt(arena_A.val());
  arena_t<ret_type> res = res_val;

//   arena_t<promote_scalar_t<double, T>> res_val
//       = Eigen::MatrixXd::Identity(A.matrix().rows(), A.matrix().cols());

//  arena_A.ldlt().solveInPlace(res_val);

//  arena_t<ret_type> res = res_val;

  reverse_pass_callback([res, res_val, arena_A]() mutable {
    arena_A.adj() -= (res_val * res.val().eval())
                           .template triangularView<Eigen::Lower>();
   // arena_A.adj() -= res_val.transpose() * res.adj_op() * res_val.transpose();
  });

  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif
