#ifndef STAN_MATH_PRIM_FUN_MEAN_HPP
#define STAN_MATH_PRIM_FUN_MEAN_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the sample mean (i.e., average) of the coefficients
 * in the specified standard vector.
 *
 * @tparam T type of elements in the vector
 * @param v Specified vector.
 * @return Sample mean of vector coefficients.
 * @throws std::domain_error if the size of the vector is less
 * than 1.
 */
template <typename StdVec, require_std_vector_t<StdVec>* = nullptr,
          require_not_vt_var<StdVec>* = nullptr>
inline value_type_t<StdVec> mean(StdVec&& v) {
  using vec_value = value_type_t<StdVec>;
  check_nonzero_size("mean", "v", v);
  Eigen::Map<const Eigen::Matrix<vec_value, Eigen::Dynamic, 1>> m(&v[0],
                                                                  v.size());
  return m.mean();
}

/**
 * Returns the sample mean (i.e., average) of the coefficients
 * in the specified vector, row vector, or matrix.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m Specified vector, row vector, or matrix.
 * @return Sample mean of vector coefficients.
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_vt_var<EigMat>* = nullptr>
inline value_type_t<EigMat> mean(EigMat&& m) {
  check_nonzero_size("mean", "m", m);
  return m.mean();
}

}  // namespace math
}  // namespace stan

#endif
