#ifndef STAN_MATH_PRIM_FUN_VARIANCE_HPP
#define STAN_MATH_PRIM_FUN_VARIANCE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/mean.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the sample variance (divide by length - 1) of the
 * coefficients in the specified standard vector.
 *
 * @tparam T type of elements in the vector
 * @param v specified vector
 * @return sample variance of vector
 * @throw <code>std::invalid_argument</code> if the vector has size zero
 */
template <typename StdVec, require_std_vector_t<StdVec>* = nullptr,
          require_not_vt_var<StdVec>* = nullptr>
inline value_type_t<StdVec> variance(StdVec&& v) {
  using vec_value = value_type_t<StdVec>;
  check_nonzero_size("variance", "v", v);
  if (v.size() == 1) {
    return vec_value{0.0};
  }
  const vec_value v_mean(mean(v));
  vec_value sum_sq_diff(0);
  for (size_t i = 0; i < v.size(); ++i) {
    const vec_value diff = v[i] - v_mean;
    sum_sq_diff += diff * diff;
  }
  return sum_sq_diff / (v.size() - 1);
}

/**
 * Returns the sample variance (divide by length - 1) of the
 * coefficients in the specified matrix
 *
 * @tparam T type of elements in the vector
 * @tparam R number of rows in the matrix, can be Eigen::Dynamic
 * @tparam C number of columns in the matrix, can be Eigen::Dynamic
 *
 * @param m matrix
 * @return sample variance of coefficients
 * @throw <code>std::invalid_argument</code> if the matrix has size zero
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_vt_var<EigMat>* = nullptr>
inline value_type_t<EigMat> variance(EigMat&& m) {
  using eig_value = value_type_t<EigMat>;
  using ref_inner = const typename std::decay_t<EigMat>::PlainObject;
  check_nonzero_size("variance", "m", m);
  if (m.size() == 1) {
    return eig_value{0.0};
  }
  const Eigen::Ref<ref_inner>& mat = m;
  return (mat.array() - mat.mean()).square().sum()
         / eig_value(mat.size() - 1.0);
}

}  // namespace math
}  // namespace stan

#endif
