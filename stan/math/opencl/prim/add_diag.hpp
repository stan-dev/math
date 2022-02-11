#ifndef STAN_MATH_OPENCL_PRIM_ADD_DIAG_HPP
#define STAN_MATH_OPENCL_PRIM_ADD_DIAG_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Returns a Matrix with values added along the main diagonal
 *
 * @tparam T_m type of input kernel generator expression for the input matrix
 * @tparam T_a type of input kernel generator expression to add along the
 * diagonal
 *
 * @param mat input kernel generator expression
 * @param to_add scalar value or input kernel generator expression to add along
 * the diagonal
 * @return a kernel generator expressio with to_add added along main diagonal
 */
template <typename T_m, typename T_a,
          require_all_kernel_expressions_and_none_scalar_t<T_m>* = nullptr,
          require_all_kernel_expressions_t<T_a>* = nullptr>
inline auto add_diag(T_m&& mat, T_a&& to_add) {  // NOLINT
  if (!is_stan_scalar<T_a>::value) {
    const size_t length_diag = std::min(mat.rows(), mat.cols());
    check_consistent_sizes("add_diag (OpenCL)", "number of elements of to_add",
                           to_add, "diagonal", length_diag);
  }
  matrix_cl<value_type_t<T_m>> mat_eval = mat;
  diagonal(mat_eval) = diagonal(mat_eval) + to_add;
  return mat_eval;
}
}  // namespace math
}  // namespace stan

#endif
#endif
