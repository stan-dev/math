#ifndef STAN_MATH_PRIM_FUN_SIZE_MVT_HPP
#define STAN_MATH_PRIM_FUN_SIZE_MVT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stdexcept>
#include <vector>

namespace stan {
namespace math {

/** \ingroup type_trait
 * Provides the size of a multivariate argument.
 *
 * This is the default template function. For any scalar type, this
 * will throw an std::invalid_argument exception since a scalar is not
 * a multivariate structure.
 *
 * @tparam T type to take size of. The default template function should
 *   only match scalars.
 * @throw std::invalid_argument since the type is a scalar.
 */
template <typename ScalarT, require_stan_scalar_t<ScalarT>* = nullptr>
size_t size_mvt(const ScalarT& /* unused */) {
  throw std::invalid_argument("size_mvt passed to an unrecognized type.");
}

template <typename MatrixT, require_matrix_t<MatrixT>* = nullptr>
size_t size_mvt(const MatrixT& /* unused */) {
  return 1U;
}

template <typename MatrixT, require_matrix_t<MatrixT>* = nullptr>
size_t size_mvt(const std::vector<MatrixT>& x) {
  return x.size();
}

}  // namespace math
}  // namespace stan
#endif
