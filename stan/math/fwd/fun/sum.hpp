#ifndef STAN_MATH_FWD_FUN_SUM_HPP
#define STAN_MATH_FWD_FUN_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/fwd/core.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the sum of the entries of the specified standard
 * vector.
 *
 * @tparam T type of elements in the vector
 * @param m Vector.
 * @return Sum of vector entries.
 */
template <typename T>
inline fvar<T> sum(const std::vector<fvar<T>>& m) {
  if (m.size() == 0) {
    return 0.0;
  }
  std::vector<T> vals(m.size());
  std::vector<T> tans(m.size());
  for (size_t i = 0; i < m.size(); ++i) {
    vals[i] = m[i].val();
    tans[i] = m[i].tangent();
  }
  return fvar<T>(sum(vals), sum(tans));
}

/**
 * Return the sum of the entries of the specified matrix.
 *
 * @tparam T type of the matrix
 *
 * @param m Matrix.
 * @return Sum of matrix entries.
 */
template <typename T, require_eigen_vt<is_fvar, T>* = nullptr>
inline value_type_t<T> sum(const T& m) {
  if (m.size() == 0) {
    return 0.0;
  }
  const Eigen::Ref<const plain_type_t<T>>& m_ref = m;
  return {sum(m_ref.val()), sum(m_ref.d())};
}

}  // namespace math
}  // namespace stan
#endif
