#ifndef STAN_MATH_REV_FUN_DOT_SELF_HPP
#define STAN_MATH_REV_FUN_DOT_SELF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the dot product of a vector of var with itself.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase and
 * have one compile time dimension equal to 1)
 * @param[in] v Vector.
 * @return Dot product of the vector with itself.
 */
template <typename T, require_eigen_vector_vt<is_var, T>* = nullptr>
inline var dot_self(const T& v) {
  arena_matrix<plain_type_t<decltype(value_of(v))>> v_val = value_of(v);
  arena_matrix<plain_type_t<T>> arena_v = v;

  var res = v_val.dot(v_val);

  reverse_pass_callback(
      [=]() mutable { arena_v.adj() += 2.0 * res.adj() * v_val; });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
