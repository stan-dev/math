#ifndef STAN_MATH_REV_FUN_VARIANCE_HPP
#define STAN_MATH_REV_FUN_VARIANCE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the sample variance of the specified standard
 * vector, Eigen vector, Eigen row vector, or Eigen matrix
 * Raise domain error if size is not greater than zero.
 *
 * @param[in] x a input
 * @return sample variance of input
 */
template <typename T, require_container_vt<is_var, T>* = nullptr>
inline var variance(const T& x) {
  check_nonzero_size("variance", "x", x);
  if (x.size() == 1) {
    return 0;
  }

  const auto& x_ref = to_ref(x);
  const auto& x_array = as_array_or_scalar(x_ref);
  arena_t<decltype(x_array)> arena_x = x_array;
  
  const auto& x_val = to_ref(value_of(x_array));
  double mean = x_val.mean();
  arena_matrix<plain_type_t<decltype(x_val)>>
    arena_diff = x_val.array() - mean;

  var res = arena_diff.matrix().squaredNorm() / (x_array.size() - 1);

  reverse_pass_callback([arena_x, res, arena_diff]() mutable {
    arena_x.adj() += (2.0 * res.adj() / (arena_x.size() - 1)) *
      arena_diff;
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
