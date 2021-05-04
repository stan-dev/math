#ifndef STAN_MATH_PRIM_FUN_UNIT_VECTOR_FREE_HPP
#define STAN_MATH_PRIM_FUN_UNIT_VECTOR_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Transformation of a unit length vector to a "free" vector
 * However, we are just fixing the unidentified radius to 1.
 * Thus, the transformation is just the identity
 *
 * @tparam EigVec A type derived from `EigenBase` with 1 compile time row or
 * column.
 * @param x unit vector of dimension K
 * @return Unit vector of dimension K considered "free"
 */
template <typename EigVec, require_eigen_col_vector_t<EigVec>* = nullptr>
inline auto unit_vector_free(EigVec&& x) {
  auto&& x_ref = to_ref(std::forward<EigVec>(x));
  check_unit_vector("stan::math::unit_vector_free", "Unit vector variable",
                    x_ref);
  return x_ref;
}

}  // namespace math
}  // namespace stan

#endif
