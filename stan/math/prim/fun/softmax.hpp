#ifndef STAN_MATH_PRIM_FUN_SOFTMAX_HPP
#define STAN_MATH_PRIM_FUN_SOFTMAX_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the softmax of the specified vector.
 *
 * <p>
 * \f$
 * \mbox{softmax}(y)
 * = \frac{\exp(y)}
 * {\sum_{k=1}^K \exp(y_k)},
 * \f$
 *
 * <p>The entries in the Jacobian of the softmax function are given by
 * \f$
 * \begin{array}{l}
 * \displaystyle
 * \frac{\partial}{\partial y_m} \mbox{softmax}(y)[k]
 * \\[8pt]
 * \displaystyle
 * \mbox{ } \ \ \ = \left\{
 * \begin{array}{ll}
 * \mbox{softmax}(y)[k] \times (1 - \mbox{softmax}(y)[m])
 * & \mbox{ if } m = k, \mbox{ and}
 * \\[6pt]
 * -\mbox{softmax}(y)[k] \times \mbox{softmax}(y)[m]
 * & \mbox{ if } m \neq k.
 * \end{array}
 * \right.
 * \end{array}
 * \f$
 *
 * @tparam Vec type of elements in the vector
 * @param[in] v Vector to transform.
 * @return Unit simplex result of the softmax transform of the vector.
 */
template <typename Vec,
          require_eigen_vector_vt<std::is_arithmetic, Vec>* = nullptr>
inline plain_type_t<Vec> softmax(const Vec& v) {
  if (v.size() == 0) {
    return v;
  }
  const auto& v_ref = to_ref(v);
  const auto theta = (v_ref.array() - v_ref.maxCoeff()).exp().eval();
  return (theta / theta.sum()).matrix();
}

/**
 * Return the softmax of the rows of the specified matrix.
 * Each row is transformed independently; the result is a row-stochastic
 * matrix whose rows each sum to one.
 *
 * @tparam Mat type of input matrix
 * @param[in] m Matrix to transform row-wise.
 * @return Row-stochastic matrix result of applying softmax to each row.
 */
template <typename Mat, require_eigen_vt<std::is_arithmetic, Mat>* = nullptr,
          require_not_eigen_vector_t<Mat>* = nullptr>
inline plain_type_t<Mat> softmax(const Mat& m) {
  const auto& m_ref = to_ref(m);
  const auto shifted
      = (m_ref.array().colwise() - m_ref.rowwise().maxCoeff().array()).eval();
  const auto exp_s = shifted.exp().eval();
  return (exp_s.colwise() / exp_s.rowwise().sum()).matrix();
}

}  // namespace math
}  // namespace stan

#endif
