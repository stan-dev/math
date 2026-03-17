#ifndef STAN_MATH_REV_FUN_SUM_HPP
#define STAN_MATH_REV_FUN_SUM_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the sum of the entries of the specified vector.
 *
 * @param m Vector.
 * @return Sum of vector entries.
 */
template <typename Alloc>
inline var sum(const std::vector<var, Alloc>& m) {
  if (unlikely(m.empty())) {
    return 0.0;
  } else {
    auto arena_m = to_arena(as_array_or_scalar(m));
    return make_callback_var(arena_m.val().sum(), [arena_m](auto& vi) mutable {
      arena_m.adj() += vi.adj();
    });
  }
}

/**
 * Returns the sum of the coefficients of the specified
 * matrix.
 *
 * @tparam T type of the matrix of vector. Can be either a var matrix or
 *  matrix of vars.
 * @param x Specified var_value containing a matrix or vector.
 * @return Sum of coefficients of matrix.
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline var sum(T&& x) {
  arena_t<T> x_arena(std::forward<T>(x));
  return make_callback_var(sum(x_arena.val()), [x_arena](auto& vi) mutable {
    x_arena.adj().array() += vi.adj();
  });
}

}  // namespace math
}  // namespace stan
#endif
