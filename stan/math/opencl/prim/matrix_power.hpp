#ifndef STAN_MATH_OPENCL_PRIM_MATRIX_POWER_HPP
#define STAN_MATH_OPENCL_PRIM_MATRIX_POWER_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/diag_matrix.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Returns the nth power of the specific matrix. M^n = M * M * ... * M.
 *
 * @tparam T type of the matrix
 *
 * @param[in] M a square matrix
 * @param[in] n exponent
 * @return nth power of M
 * @throw std::domain_error if the matrix contains NaNs or infinities.
 * @throw std::invalid_argument if the exponent is negative or the matrix is not
 * square.
 */
template <typename T_m,
          require_all_kernel_expressions_and_none_scalar_t<T_m>* = nullptr>
inline plain_type_t<T_m> matrix_power(T_m&& M, const int n) {
  const char* function = "matrix_power(OpenCL)";
  check_square(function, "M", M);
  check_nonnegative(function, "n", n);
  plain_type_t<T_m> MM = std::forward<T_m>(M);
  check_cl(function, "M", value_of(MM), "finite") = isfinite(value_of(MM));
  if (n == 0)
    return diag_matrix(constant(1.0, M.rows(), 1));
  plain_type_t<T_m> result = MM;
  for (int nn = n - 1; nn > 0; nn /= 2) {
    if (nn % 2 == 1) {
      result = result * MM;
      --nn;
    }
    MM = MM * MM;
  }
  return result;
}
}  // namespace math
}  // namespace stan

#endif
#endif
