#ifndef STAN_MATH_PRIM_FUN_SINGULAR_VALUES_HPP
#define STAN_MATH_PRIM_FUN_SINGULAR_VALUES_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>

namespace stan {
namespace math {

/**
 * Return the vector of the singular values of the specified matrix
 * in decreasing order of magnitude.
 * <p>See the documentation for <code>svd()</code> for
 * information on the singular values.
 *
 * @tparam EigMat type of the matrix
 * @param m Specified matrix.
 * @return Singular values of the matrix.
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr,
          require_not_st_var<EigMat>* = nullptr>
Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, 1> singular_values(
    const EigMat& m) {
  check_nonzero_size("singular_values", "m", m);

  return Eigen::JacobiSVD<Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic,
                                        Eigen::Dynamic> >(m)
      .singularValues();
}

}  // namespace math
}  // namespace stan

#endif
