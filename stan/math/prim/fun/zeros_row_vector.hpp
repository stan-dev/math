#ifndef STAN_MATH_PRIM_FUN_ZEROS_ROW_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_ZEROS_ROW_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a row vector of zeros
 *
 * @param K size of the row vector
 * @return A row vector of size K with all elements initialized to 0.
 * @throw std::domain_error if K is negative.
 */
inline auto zeros_row_vector(int K) {
  check_nonnegative("zeros_row_vector", "size", K);
  return Eigen::RowVectorXd::Zero(K);
}

}  // namespace math
}  // namespace stan

#endif
