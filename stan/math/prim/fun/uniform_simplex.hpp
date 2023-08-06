#ifndef STAN_MATH_PRIM_FUN_UNIFORM_SIMPLEX_HPP
#define STAN_MATH_PRIM_FUN_UNIFORM_SIMPLEX_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a uniform simplex of size K
 *
 * @param K size of the simplex
 * @return A vector of size K with all elements initialized to 1 / K,
 * so that their sum is equal to 1.
 * @throw std::domain_error if K is not positive.
 */
inline auto uniform_simplex(int K) {
  check_positive("uniform_simplex", "size", K);
  return Eigen::VectorXd::Constant(K, 1.0 / K);
}

}  // namespace math
}  // namespace stan

#endif
