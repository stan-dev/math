#ifndef STAN_MATH_OPENCL_REV_TO_MATRIX_HPP
#define STAN_MATH_OPENCL_REV_TO_MATRIX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Returns a matrix representation of a vector or matrix in column-major
 * order with the specified number of rows and columns.
 *
 * @tparam T_x type of the matrix
 *
 * @param x matrix
 * @param m rows
 * @param n columns
 * @return Reshaped input matrix
 * @throw <code>std::invalid_argument</code> if the sizes
 * do not match
 */
template <typename T_x,
          require_all_kernel_expressions_and_none_scalar_t<T_x>* = nullptr>
inline var_value<matrix_cl<double>> to_matrix(const var_value<T_x>& x, int m,
                                              int n) {
  return make_callback_var(
      to_matrix(value_of(x), m, n),
      [x, m, n](vari_value<matrix_cl<double>>& res) mutable {
        matrix_cl<double> x_adj_cpy = std::move(x.adj());
        matrix_cl<double> reshaped(x_adj_cpy.buffer(), m, n);
        for (cl::Event e : x_adj_cpy.read_events()) {
          reshaped.add_read_event(e);
        }
        for (cl::Event e : x_adj_cpy.write_events()) {
          reshaped.add_write_event(e);
        }
        reshaped += res.adj();
        for (cl::Event e : reshaped.write_events()) {
          x_adj_cpy.add_write_event(e);
        }
        x.adj() = std::move(x_adj_cpy);
      });
}

}  // namespace math
}  // namespace stan
#endif
#endif
