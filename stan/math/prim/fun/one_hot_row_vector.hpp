#ifndef STAN_MATH_PRIM_FUN_ONE_HOT_ROW_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_ONE_HOT_ROW_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a row vector with 1 in the k-th position and zero elsewhere
 *
 * @param K size of the row vector
 * @param k position of the 1 (indexing from 1)
 * @return A row vector of size K with all elements initialized to zero
 * and a 1 in the k-th position.
 * @throw std::domain_error if K is not positive, or if k is less than 1 or
 * greater than K.
 */
inline Eigen::RowVectorXd one_hot_row_vector(int K, int k) {
  static const char* function = "one_hot_row_vector";
  check_positive(function, "size", K);
  check_bounded(function, "k", k, 1, K);

  Eigen::RowVectorXd ret = Eigen::RowVectorXd::Zero(K);
  ret(k - 1) = 1;
  return ret;
}

}  // namespace math
}  // namespace stan

#endif
