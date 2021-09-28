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
 * @tparam EigVec type of the vector (must be derived from \c Eigen::MatrixBase
 * and have one compile time dimension equal to 1)
 *
 * @param m Vector of values.
 * @return Cumulative sum of values.
 */
template <typename EigVec, require_rev_matrix_t<EigVec>* = nullptr>
inline auto cumulative_sum(const EigVec& m) {
  if (m.size() < 2) {
    return m;
  }
  arena_t<EigVec> m_arena = m;
  arena_t<EigVec> ret = cumulative_sum(m_arena.val());
  reverse_pass_callback([m_arena, ret]() mutable {
    for (Eigen::Index i = ret.size() - 1; i > 0; --i) {
      ret.adj()(i - 1) += ret.adj()(i);
    }
    m_arena.adj() += ret.adj();
  });
  return plain_type_t<EigVec>(ret);
}

}  // namespace math
}  // namespace stan

#endif
