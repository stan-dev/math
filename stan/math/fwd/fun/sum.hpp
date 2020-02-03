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
 * @param m Vector.
 * @return Sum of vector entries.
 */
template <typename T>
inline fvar<T> sum(const std::vector<fvar<T> >& m) {
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
 * @tparam T inner type of the fvar matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m Matrix.
 * @return Sum of matrix entries.
 */
template <typename T, int R, int C>
inline fvar<T> sum(const Eigen::Matrix<fvar<T>, R, C>& m) {
  if (m.size() == 0) {
    return 0.0;
  }
  Eigen::Matrix<T, Eigen::Dynamic, 1> vals(m.size());
  Eigen::Matrix<T, Eigen::Dynamic, 1> tans(m.size());
  for (int i = 0; i < m.size(); ++i) {
    vals(i) = m(i).val();
    tans(i) = m(i).tangent();
  }
  return fvar<T>(sum(vals), sum(tans));
}

}  // namespace math
}  // namespace stan
#endif
