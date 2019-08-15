#ifndef STAN_MATH_PRIM_MAT_FUN_VARIANCE_HPP
#define STAN_MATH_PRIM_MAT_FUN_VARIANCE_HPP

#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/fun/mean.hpp>
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
inline auto variance(const std::vector<T>& v) {
  check_nonzero_size("variance", "v", v);
  if (v.size() == 1)
    return scalar_type_t<T>(0.0);
  T v_mean(mean(v));
  T sum_sq_diff(0);
  for (size_t i = 0; i < v.size(); ++i) {
    T diff = v[i] - v_mean;
    sum_sq_diff += diff * diff;
  }
  return sum_sq_diff / (v.size() - 1);
}

/**
 * Returns the sample variance (divide by length - 1) of the
 * coefficients in the specified column vector.
 * @param m Specified vector.
 * @return Sample variance of vector.
 */
template <typename T, typename = enable_if_eigen<T>>
inline auto variance(const T& m) {
  check_nonzero_size("variance", "m", m);

  if (m.size() == 1)
    return typename T::Scalar(0.0);
  return (m.array() - mean(m)).square().sum() / (m.size() - 1);
}

}  // namespace math
}  // namespace stan
#endif
