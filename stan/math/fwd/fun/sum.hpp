#ifndef STAN_MATH_FWD_FUN_SUM_HPP
#define STAN_MATH_FWD_FUN_SUM_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the sum of the entries of the specified standard
 * vector.
 *
 * @tparam T type of elements in the vector
 * @param v Vector.
 * @return Sum of vector entries.
 */
template <typename T>
inline fvar<T> sum(const std::vector<fvar<T>>& v) {
  if (v.size() == 0) {
    return 0.0;
  }
  std::vector<T> vals(v.size());
  std::vector<T> tans(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    vals[i] = v[i].val();
    tans[i] = v[i].tangent();
  }
  return fvar<T>(sum(vals), sum(tans));
}

/**
 * Return the sum of the entries of the specified matrix.
 *
 * @tparam T type of the matrix
 *
 * @param v Matrix.
 * @return Sum of matrix entries.
 */
template <typename T, require_eigen_vt<is_fvar, T>* = nullptr>
inline value_type_t<T> sum(const T& v) {
  if (v.size() == 0) {
    return 0.0;
  }
  const Eigen::Ref<const plain_type_t<T>>& v_ref = v;
  return {sum(v_ref.val()), sum(v_ref.d())};
}

}  // namespace math
}  // namespace stan
#endif
