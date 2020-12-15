#ifndef STAN_MATH_OPENCL_REV_COL_HPP
#define STAN_MATH_OPENCL_REV_COL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Return the specified column of the specified 
 * `var_value<matrix_cl<double>>` using start-at-1 indexing.
 *
 * @param m Matrix.
 * @param j Column index (count from 1).
 * @return Specified column of the matrix.
 * @throw std::out_of_range if j is out of range.
 */
inline var_value<matrix_cl<double>> col(
    const var_value<matrix_cl<double>>& m, size_t j) {
  var_value<matrix_cl<double>> res = col(m.val(), j);

  reverse_pass_callback([m, res, j]() mutable {
    block(m.adj(), 0, j-1, m.rows(), 1) = block(m.adj(), 0, j-1, m.rows(), 1) + res.adj();
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
