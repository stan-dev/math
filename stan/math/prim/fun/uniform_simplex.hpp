#ifndef STAN_MATH_PRIM_FUN_UNIFORM_SIMPLEX_HPP
#define STAN_MATH_PRIM_FUN_UNIFORM_SIMPLEX_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a uniform simplex of length K
 *
 * @param K length of the simplex
 * @return A vector of length K with all elements initialised to 1 / K.
 */
Eigen::VectorXd uniform_simplex(int K) {
  check_nonnegative("uniform_simplex", "length", K);
  return Eigen::VectorXd::Constant(K, 1.0 / K);
}

}  // namespace math
}  // namespace stan

#endif
