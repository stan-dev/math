#ifndef STAN_MATH_PRIM_FUN_ONES_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_ONES_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a vector of ones
 *
 * @param K length of the vector
 * @return A vector of length K with all elements initialised to 1.
 */
Eigen::VectorXd ones_vector(int K) {
  check_nonnegative("ones_vector", "length", K);
  return Eigen::VectorXd::Constant(K, 1);
}

}  // namespace math
}  // namespace stan

#endif
