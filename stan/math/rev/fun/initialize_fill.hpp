#ifndef STAN_MATH_REV_FUN_INITIALIZE_FILL_HPP
#define STAN_MATH_REV_FUN_INITIALIZE_FILL_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/initialize_fill.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Fill the specified container with the specified value. This function does
 * not perform a callback to propogate the adjoints upward
 *
 * The specified matrix is filled by element.
 *
 * @tparam VarMat a `var_value` with inner type from `EigenBase`
 * @tparam S A var.
 *
 * @param x Container.
 * @param y Value.
 */
template <typename VarMat, typename S, require_var_matrix_t<VarMat>* = nullptr,
          require_stan_scalar_t<S>* = nullptr>
inline void initialize_fill(VarMat& x, const S& y) {
  x.vi_->val_.fill(value_of(y));
}

}  // namespace math
}  // namespace stan

#endif
