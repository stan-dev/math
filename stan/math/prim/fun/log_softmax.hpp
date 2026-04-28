#ifndef STAN_MATH_PRIM_FUN_LOG_SOFTMAX_HPP
#define STAN_MATH_PRIM_FUN_LOG_SOFTMAX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the softmax of the specified
 * vector.
 *
 * \f$
 * \log \mbox{softmax}(y)
 * \ = \ y - \log \sum_{k=1}^K \exp(y_k)
 * \ = \ y - \mbox{log\_sum\_exp}(y).
 * \f$
 *
 * For the log softmax function, the entries in the Jacobian are
 * \f$
 * \frac{\partial}{\partial y_m} \mbox{softmax}(y)[k]
 * = \left\{
 * \begin{array}{ll}
 * 1 - \mbox{softmax}(y)[m]
 * & \mbox{ if } m = k, \mbox{ and}
 * \\[6pt]
 * \mbox{softmax}(y)[m]
 * & \mbox{ if } m \neq k.
 * \end{array}
 * \right.
 * \f$
 *
 * @tparam Container type of input vector to transform
 * @param[in] x vector to transform
 * @return log unit simplex result of the softmax transform of the vector.
 */
template <typename Container, require_st_arithmetic<Container>* = nullptr,
          require_container_t<Container>* = nullptr,
          require_not_t<bool_constant<
              is_eigen<std::decay_t<Container>>::value
              && !is_eigen_vector<std::decay_t<Container>>::value>>* = nullptr>
inline auto log_softmax(Container&& x) {
  check_nonzero_size("log_softmax", "v", x);
  return make_holder(
      [](auto&& a) {
        return apply_vector_unary<ref_type_t<Container>>::apply(
            std::forward<decltype(a)>(a),
            [](auto&& v) { return v.array() - log_sum_exp(v); });
      },
      to_ref(std::forward<Container>(x)));
}

/**
 * Return the log softmax of the rows of the specified matrix.
 * Each row is transformed independently; the result has the same shape
 * as the input.
 *
 * @tparam Mat type of input matrix
 * @param[in] m Matrix to transform row-wise.
 * @return Log-softmax applied row-wise.
 */
template <typename Mat, require_eigen_vt<std::is_arithmetic, Mat>* = nullptr,
          require_not_eigen_vector_t<Mat>* = nullptr>
inline plain_type_t<Mat> log_softmax(const Mat& m) {
  check_nonzero_size("log_softmax", "m", m);
  const auto& m_ref = to_ref(m);
  const auto shifted
      = (m_ref.array().colwise() - m_ref.rowwise().maxCoeff().array()).eval();
  const auto exp_s = shifted.exp().eval();
  return (shifted.colwise() - exp_s.rowwise().sum().log()).matrix();
}

}  // namespace math
}  // namespace stan
#endif
