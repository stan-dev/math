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
 * coefficients in the specified matrix
 *
 * @tparam EigMat type inheriting from `EigenBase` that does not have an `var`
 *  scalar type.
 *
 * @param m matrix
 * @return sample variance of coefficients
 * @throw <code>std::invalid_argument</code> if the matrix has size zero
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_vt_var<EigMat>* = nullptr>
inline value_type_t<EigMat> variance(const EigMat& m) {
  using value_t = value_type_t<EigMat>;
  const auto& mat = to_ref(m);
  check_nonzero_size("variance", "m", mat);
  if (mat.size() == 1) {
    return value_t{0.0};
  }
  return (mat.array() - mat.mean()).square().sum() / value_t(mat.size() - 1.0);
}

/**
 * Returns the sample variance (divide by length - 1) of the
 * coefficients in the specified standard vector.
 *
 * @tparam StdVec A standard library vector that does not contain a var.
 * @param v specified vector
 * @return sample variance of vector
 * @throw <code>std::invalid_argument</code> if the vector has size zero
 */
template <typename StdVec, require_std_vector_t<StdVec>* = nullptr,
          require_not_vt_var<StdVec>* = nullptr>
inline value_type_t<StdVec> variance(const StdVec& v) {
  using eigen_t = Eigen::Matrix<value_type_t<StdVec>, -1, 1>;
  return variance(Eigen::Map<const eigen_t>(v.data(), v.size()));
}

}  // namespace math
}  // namespace stan

#endif
