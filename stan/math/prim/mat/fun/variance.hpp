#ifndef STAN_MATH_PRIM_MAT_FUN_VARIANCE_HPP
#define STAN_MATH_PRIM_MAT_FUN_VARIANCE_HPP

#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <boost/math/tools/promotion.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the sample variance (divide by length - 1) of the
 * coefficients in the specified standard vector.
 * @param v Specified vector.
 * @return Sample variance of vector.
 * @throws std::domain_error if the size of the vector is less
 * than 1.
 */
template <typename T>
inline typename boost::math::tools::promote_args<T>::type variance(
    const std::vector<T>& v) {
  check_nonzero_size("variance", "v", v);
  if (v.size() == 1)
    return 0.0;
  Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> m(&v[0], v.size());

  return (m.array() - m.mean()).matrix().squaredNorm() / (m.size() - 1);
}

/**
 * Returns the sample variance (divide by length - 1) of the
 * coefficients in the specified column vector.
 * @param m Specified vector.
 * @return Sample variance of vector.
 */
template <typename T, int R, int C>
inline typename boost::math::tools::promote_args<T>::type variance(
    const Eigen::Matrix<T, R, C>& m) {
  check_nonzero_size("variance", "m", m);

  if (m.size() == 1)
    return 0.0;

  return (m.array() - m.mean()).matrix().squaredNorm() / (m.size() - 1);
}

}  // namespace math
}  // namespace stan
#endif
