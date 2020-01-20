#ifndef STAN_MATH_PRIM_FUN_ELT_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_ELT_MULTIPLY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise multiplication of the specified
 * matrices.
 *
 * @tparam T1 type of the first matrix or expression
 * @tparam T2 type of the second matrix or expression
 *
 * @param m1 First matrix or expression
 * @param m2 Second matrix or expression
 * @return Elementwise product of matrices.
 */
template <typename T1, typename T2, typename = require_all_eigen_t<T1, T2>>
auto elt_multiply(const T1& m1, const T2& m2) {
  check_matching_dims("elt_multiply", "m1", m1, "m2", m2);
  return m1.cwiseProduct(m2);
}

}  // namespace math
}  // namespace stan

#endif
