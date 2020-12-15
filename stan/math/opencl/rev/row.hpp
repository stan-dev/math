#ifndef STAN_MATH_OPENCL_REV_ROW_HPP
#define STAN_MATH_OPENCL_REV_ROW_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Return the specified row of the specified 
 * `var_value<matrix_cl<double>>` using start-at-1 indexing.
 *
 * @param m Matrix.
 * @param j Row index (count from 1).
 * @return Specified row of the matrix.
 * @throw std::out_of_range if j is out of range.
 */
inline var_value<matrix_cl<double>> row(
    const var_value<matrix_cl<double>>& m, size_t j) {
  var_value<matrix_cl<double>> res = row(m.val(), j);

  reverse_pass_callback([m, res, j]() mutable {
    block(m.adj(), j-1, 0, 1, m.rows()) = block(m.adj(), j-1, 0, 1, m.rows()) + res.adj();
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
