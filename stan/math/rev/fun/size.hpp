#ifndef STAN_MATH_REV_FUN_SIZE_HPP
#define STAN_MATH_REV_FUN_SIZE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core/var.hpp>
#include <cstdlib>
#include <vector>

namespace stan {
namespace math {

/** \ingroup type_trait
 * Returns the length of var_value<double> (always of length 1)
 */
inline size_t size(const var_value<double>& /*x*/) {
  return 1U;
}

/** \ingroup type_trait
 * Returns the size of the provided var_value<Eigen::Matrix<T, R, C>>
 *
 * @param m var_value
 * @tparam T type of m
 */
template<int R, int C>
inline size_t size(const var_value<Eigen::Matrix<double, R, C>>& m) {
  return m.val().size();
}

}  // namespace math
}  // namespace stan
#endif
