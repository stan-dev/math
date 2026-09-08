#ifndef STAN_MATH_REV_FUN_INVERSE_HPP
#define STAN_MATH_REV_FUN_INVERSE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/inverse.hpp>

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

  using ret_type = return_var_matrix_t<T>;
  if (unlikely(m.size() == 0)) {
    return ret_type(m);
  }

  arena_t<T> arena_m = m;
  arena_t<promote_scalar_t<double, T>> res_val = arena_m.val().inverse();
  arena_t<ret_type> res = res_val;

  reverse_pass_callback([res, res_val, arena_m]() mutable {
    arena_m.adj() -= res_val.transpose() * res.adj_op() * res_val.transpose();
  });

  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif
