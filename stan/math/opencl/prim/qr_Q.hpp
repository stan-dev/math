#ifndef STAN_MATH_OPENCL_PRIM_QR_Q_HPP
#define STAN_MATH_OPENCL_PRIM_QR_Q_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/qr_decomposition.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Returns the orthogonal factor of the fat QR decomposition
 *
 * @tparam T_m type of the matrix
 * @param m Matrix.
 * @return Orthogonal matrix with maximal columns
 */
template <typename T_m,
          require_all_kernel_expressions_and_none_scalar_t<T_m>* = nullptr>
inline matrix_cl<double> qr_Q(T_m&& m) {
  check_nonzero_size("qr_Q(OpenCL)", "m", m);

  matrix_cl<double> mat_eval = std::forward<T_m>(m);
  matrix_cl<double> Q, R;

  qr_decomposition_cl(m, Q, R);

  auto Q_left
      = block_zero_based(Q, 0, 0, Q.rows(), std::min(m.rows(), m.cols()));
  Q_left = select(colwise_broadcast(transpose(diagonal(R) < 0.0)), -Q_left,
                  Q_left);
  return Q;
}
}  // namespace math
}  // namespace stan

#endif
#endif
