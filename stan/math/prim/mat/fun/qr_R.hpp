#ifndef STAN_MATH_PRIM_MAT_FUN_QR_R_HPP
#define STAN_MATH_PRIM_MAT_FUN_QR_R_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Returns the upper triangular factor of the fat QR decomposition
 * @param m Matrix.
 * @tparam T scalar type
 * @return Upper triangular matrix with maximal rows
 */
template <typename T, typename = enable_if_eigen<T>>
auto qr_R(const T& m) {
  check_nonzero_size("qr_R", "m", m);
  check_greater_or_equal("qr_R", "m.rows()", m.rows(), m.cols());

  using plain_type = typename T::PlainObject;
  Eigen::HouseholderQR<plain_type> qr(m.rows(), m.cols());
  qr.compute(m);
  plain_type R = qr.matrixQR();
  if (m.rows() > m.cols())
    R.bottomRows(m.rows() - m.cols()).setZero();
  const int min_size = std::min(m.rows(), m.cols());
  for (int i = 0; i < min_size; i++) {
    for (int j = 0; j < i; j++)
      R.coeffRef(i, j) = 0.0;
    if (R(i, i) < 0)
      R.row(i) *= -1.0;
  }
  return R;
}
}  // namespace math
}  // namespace stan
#endif
