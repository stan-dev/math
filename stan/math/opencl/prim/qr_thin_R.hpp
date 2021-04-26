#ifndef STAN_MATH_OPENCL_PRIM_QR_THIN_R_HPP
#define STAN_MATH_OPENCL_PRIM_QR_THIN_R_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/qr_decomposition.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Returns the orthogonal factor of the thin QR decomposition
 *
 * @tparam T_m type of the matrix
 * @param m Matrix.
 * @return Orthogonal matrix with maximal columns
 */
template <typename T_m,
          require_all_kernel_expressions_and_none_scalar_t<T_m>* = nullptr>
inline matrix_cl<double> qr_thin_R(T_m&& m) {
  check_nonzero_size("qr_thin_R(OpenCL)", "m", m);

  matrix_cl<double> mat_eval = std::forward<T_m>(m);
  matrix_cl<double> Q, R;

  qr_decomposition_cl(m, Q, R);

  auto R_block
      = block_zero_based(R, 0, 0, std::min(m.rows(), m.cols()), m.cols());
  R = select(rowwise_broadcast(diagonal(R) < 0.0), -R_block, R_block).eval();
  return R;
}
}  // namespace math
}  // namespace stan

#endif
#endif
