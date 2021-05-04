#ifndef STAN_MATH_OPENCL_REV_MATRIX_POWER_HPP
#define STAN_MATH_OPENCL_REV_MATRIX_POWER_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/rev/arena_type.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/diag_matrix.hpp>
#include <stan/math/rev/core.hpp>
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
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> matrix_power(const var_value<T>& M,
                                                 const int n) {
  check_square("matrix_power", "M", M);
  check_nonnegative("matrix_power", "n", n);

  if (M.size() == 0)
    return M;

  size_t N = M.rows();
  if (n == 0) {
    return diag_matrix(constant(1.0, M.rows(), 1));
  }
  if (n == 1) {
    return M;
  }

  arena_t<std::vector<matrix_cl<double>>> arena_powers(n + 1);
  arena_powers[0] = diag_matrix(constant(1.0, M.rows(), 1));
  arena_powers[1] = M.val();
  for (size_t i = 2; i <= n; ++i) {
    arena_powers[i] = arena_powers[1] * arena_powers[i - 1];
  }

  return make_callback_var(
      arena_powers.back(),
      [M, n, arena_powers](vari_value<matrix_cl<double>> res) mutable {
        const auto& M_val = arena_powers[1];
        matrix_cl<double> adj_C = res.adj();
        matrix_cl<double> adj_M = constant(0.0, M_val.rows(), M_val.cols());
        for (size_t i = n; i > 1; --i) {
          adj_M += adj_C * transpose(arena_powers[i - 1]);
          adj_C = transpose(M_val) * adj_C;
        }
        M.adj() += adj_M + adj_C;
      });
}
}  // namespace math
}  // namespace stan

#endif
#endif
