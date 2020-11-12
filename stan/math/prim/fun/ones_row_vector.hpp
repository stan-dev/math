#ifndef STAN_MATH_PRIM_FUN_ONES_ROW_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_ONES_ROW_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a row vector of ones
 *
 * @param K size of the row vector
 * @return A row vector of size K with all elements initialised to 1.
 * @throw std::domain_error if K is negative.
 */
inline auto ones_row_vector(int K) {
  check_nonnegative("ones_row_vector", "size", K);
  return Eigen::RowVectorXd::Constant(K, 1);
}

}  // namespace math
}  // namespace stan

#endif
