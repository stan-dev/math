#ifndef STAN_MATH_PRIM_ERR_IS_COLUMN_INDEX_HPP
#define STAN_MATH_PRIM_ERR_IS_COLUMN_INDEX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if column index is in bounds.
 * By default this is a 1-indexed check (as opposed to zero-indexed).
 * Behavior can be changed by setting <code>stan::error_index::value</code>.
 * @tparam EigMat A type derived from `EigenBase`
 * @param y matrix to test
 * @param i column index to check
 * @return <code>true</code> if column index is in bounds
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline bool is_column_index(const EigMat& y, size_t i) {
  return i >= stan::error_index::value
         && i < static_cast<size_t>(y.cols()) + stan::error_index::value;
}

}  // namespace math
}  // namespace stan
#endif
