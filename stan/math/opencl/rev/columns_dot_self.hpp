#ifndef STAN_MATH_OPENCL_REV_COLUMNS_DOT_SELF_HPP
#define STAN_MATH_OPENCL_REV_COLUMNS_DOT_SELF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/size_zero.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each column of a matrix with itself.
 *
 * @tparam T type of the matrix
 * @param v Matrix.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> columns_dot_self(const var_value<T>& v) {
  if (size_zero(v)) {
    return var_value<matrix_cl<double>>(constant(0.0, 1, v.cols()));
  }

  return make_callback_var(
      columns_dot_self(v.val()),
      [v](const vari_value<matrix_cl<double>>& res) mutable {
        v.adj()
            += elt_multiply(colwise_broadcast(res.adj() * 2.0), value_of(v));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
