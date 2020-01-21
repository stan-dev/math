#ifndef STAN_MATH_PRIM_FUN_CONSTANT_ROW_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_CONSTANT_ROW_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a row vector with elements set to the same constant
 *
 * @param K size of the row vector
 * @param c constant value
 * @return A row vector of size K with all elements initialised to
 * the same constant.
 * @throw std::domain_error if K is negative.
 */
inline Eigen::RowVectorXd constant_row_vector(int K, double c) {
  check_nonnegative("constant_row_vector", "size", K);
  return Eigen::RowVectorXd::Constant(K, c);
}

}  // namespace math
}  // namespace stan

#endif
