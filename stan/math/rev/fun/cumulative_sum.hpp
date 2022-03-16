#ifndef STAN_MATH_REV_FUN_CUMULATIVE_SUM_HPP
#define STAN_MATH_REV_FUN_CUMULATIVE_SUM_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/cumulative_sum.hpp>
#include <vector>
#include <numeric>
#include <functional>

namespace stan {
namespace math {

/**
 * Return the cumulative sum of the specified vector.
 *
 * The cumulative sum is of the same type as the input and
 * has values defined by
 *
 * \code x(0), x(1) + x(2), ..., x(1) + , ..., + x(x.size()-1) @endcode
 *
 * @tparam EigVec type derived from `Eigen::EigenBase` or a `var_value<T>` with
 *  `T` deriving from `Eigen::EigenBase` with one compile time dimension
 *   equal to 1
 *
 * @param m Vector of values.
 * @return Cumulative sum of values.
 */
template <typename EigVec, require_rev_vector_t<EigVec>* = nullptr>
inline auto cumulative_sum(const EigVec& x) {
  arena_t<EigVec> x_arena(x);
  using return_t = return_var_matrix_t<EigVec>;
  arena_t<return_t> res = cumulative_sum(x_arena.val()).eval();
  if (unlikely(x.size() == 0)) {
    return return_t(res);
  }
  reverse_pass_callback([x_arena, res]() mutable {
    for (Eigen::Index i = x_arena.size() - 1; i > 0; --i) {
      x_arena.adj().coeffRef(i) += res.adj().coeffRef(i);
      res.adj().coeffRef(i - 1) += res.adj().coeffRef(i);
    }
    x_arena.adj().coeffRef(0) += res.adj().coeffRef(0);
  });
  return return_t(res);
}

}  // namespace math
}  // namespace stan

#endif
