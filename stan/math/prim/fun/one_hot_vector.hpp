#ifndef STAN_MATH_PRIM_FUN_ONE_HOT_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_ONE_HOT_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a vector with 1 in the k-th position and zero elsewhere
 *
 * @param K size of the vector
 * @param k position of the 1 (indexing from 1)
 * @return A vector of size K with all elements initialized to zero
 * and a 1 in the k-th position.
 * @throw std::domain_error if K is not positive, or if k is less than 1 or
 * greater than K.
 */
inline Eigen::VectorXd one_hot_vector(int K, int k) {
  static constexpr const char* function = "one_hot_vector";
  check_positive(function, "size", K);
  check_bounded(function, "k", k, 1, K);

  Eigen::VectorXd ret = Eigen::VectorXd::Zero(K);
  ret(k - 1) = 1;
  return ret;
}

}  // namespace math
}  // namespace stan

#endif
