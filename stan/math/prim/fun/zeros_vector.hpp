#ifndef STAN_MATH_PRIM_FUN_ZEROS_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_ZEROS_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a vector of zeros
 *
 * @param K length of the vector
 * @return A vector of length K with all elements initialised to 0.
 */
Eigen::VectorXd zeros_vector(long K) {
  check_nonnegative("ones_vector", "length", K);
  return Eigen::VectorXd::Zero(K);
}

}  // namespace math
}  // namespace stan

#endif
