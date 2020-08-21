#ifndef STAN_MATH_REV_FUN_INVERSE_HPP
#define STAN_MATH_REV_FUN_INVERSE_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/functor/arena_matrix.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>

namespace stan {
namespace math {

/**
 * Reverse mode specialization of calculating the inverse of the matrix.
 *
 * @param m specified matrix
 * @return Inverse of the matrix (an empty matrix if the specified matrix has
 * size zero).
 * @throw std::invalid_argument if the matrix is not square.
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> inverse(const T& m) {
  check_square("inverse", "m", m);
  if (m.size() == 0) {
    return {};
  }

  arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> arena_m = m;
  arena_matrix<Eigen::MatrixXd> res_val = value_of(arena_m).inverse();
  arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> res
      = res_val;

  reverse_pass_callback([res, res_val, arena_m]() mutable {
    Eigen::MatrixXd res_adj = res.adj();
    arena_m.adj() -= res_val.transpose() * res_adj * res_val.transpose();
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
