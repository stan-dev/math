#ifndef STAN_MATH_PRIM_FUN_ONES_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_ONES_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a vector of ones
 *
 * @param K size of the vector
 * @return A vector of size K with all elements initialized to 1.
 * @throw std::domain_error if K is negative.
 */
inline auto ones_vector(int K) {
  check_nonnegative("ones_vector", "size", K);
  return make_holder([](auto K_) {
   return Eigen::VectorXd::Constant(K_, 1);
  }, K);
}

}  // namespace math
}  // namespace stan

#endif
