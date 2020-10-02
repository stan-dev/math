#ifndef STAN_MATH_REV_FUN_INVERSE_HPP
#define STAN_MATH_REV_FUN_INVERSE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>

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
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto inverse(const T& m) {
  check_square("inverse", "m", m);
  if (m.size() == 0) {
    return plain_type_t<T>(m);
  }

  arena_t<T> arena_m = m;
  arena_t<T> res(value_of(arena_m).inverse().eval());
  arena_t<Eigen::MatrixXd> res_val = res.val();

  reverse_pass_callback([res, res_val, arena_m]() mutable {
    arena_m.adj().noalias() -= res_val.transpose() * res.adj_op() * res_val.transpose();
  });

  return plain_type_t<T>(res);
}

}  // namespace math
}  // namespace stan
#endif
