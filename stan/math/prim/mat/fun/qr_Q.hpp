#ifndef STAN_MATH_PRIM_MAT_FUN_QR_Q_HPP
#define STAN_MATH_PRIM_MAT_FUN_QR_Q_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Returns the orthogonal factor of the fat QR decomposition
 * @param m Matrix.
 * @tparam T scalar type
 * @return Orthogonal matrix with maximal columns
 */
template <typename T, typename = enable_if_eigen<T>>
auto qr_Q(const T& m) {
  check_nonzero_size("qr_Q", "m", m);
  check_greater_or_equal("qr_Q", "m.rows()", m.rows(), m.cols());
  using plain_type = typename T::PlainObject;
  Eigen::HouseholderQR<plain_type> qr(m.rows(), m.cols());
  qr.compute(m);
  plain_type Q = qr.householderQ();
  const int min_size = std::min(m.rows(), m.cols());
  for (int i = 0; i < min_size; i++)
    if (qr.matrixQR().coeff(i, i) < 0)
      Q.col(i) *= -1.0;
  return Q;
}

}  // namespace math
}  // namespace stan
#endif
