#ifndef STAN_MATH_REV_FUN_TCROSSPROD_HPP
#define STAN_MATH_REV_FUN_TCROSSPROD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/Eigen_NumTraits.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of post-multiplying a matrix by its
 * own transpose.
 *
 * @tparam T Type of the matrix (must be derived from \c Eigen::MatrixBase)
 * @param M Matrix to multiply.
 * @return M times its transpose.
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline plain_type_t<T> tcrossprod(const T& M) {
  using T_double = Eigen::MatrixXd;
  using T_var = plain_type_t<T>;

  if (M.size() == 0) {
    return M;
  }

  arena_t<plain_type_t<T>> arena_M = M;
  arena_t<Eigen::MatrixXd> arena_M_val = value_of(arena_M);

  arena_t<T_var> res = arena_M_val * arena_M_val.transpose();

  /*
  Eigen::MatrixXd res_val(M.rows(), M.rows());
  res_val.setZero().template selfadjointView<Eigen::Upper>().rankUpdate(arena_M_val);

  for (size_t j = 0; j < res.cols(); ++j) {
    for (size_t i = 0; i < j; ++i) {
      res.coeffRef(i, j) = res.coeffRef(j, i) = res_val.coeff(i, j);
    }
    res.coeffRef(j, j) = res_val.coeff(j, j);
    }*/

  reverse_pass_callback([res, arena_M, arena_M_val]() mutable {
    Eigen::MatrixXd adj = res.adj();
    arena_M.adj() += (adj + adj.transpose()) * arena_M_val;
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
