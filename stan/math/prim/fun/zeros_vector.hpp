#ifndef STAN_MATH_PRIM_FUN_ZEROS_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_ZEROS_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a vector of zeros.
 *
 * @param K size of the vector
 * @return A vector of size K with all elements initialised to 0.
 * @throw std::domain_error if K is negative.
 */
inline auto zeros_vector(int K) {
  check_nonnegative("zeros_vector", "size", K);
  return Eigen::VectorXd::Zero(K);
}

}  // namespace math
}  // namespace stan

#endif
