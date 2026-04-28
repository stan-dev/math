#ifndef STAN_MATH_FWD_FUN_LOG_SOFTMAX_HPP
#define STAN_MATH_FWD_FUN_LOG_SOFTMAX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun/softmax.hpp>
#include <stan/math/prim/fun/log_softmax.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Return the log softmax of the specified vector or container of vectors.
 *
 * @tparam T Type of input vector or matrix.
 * @param[in] x Unconstrained input vector.
 * @return Softmax of the input.
 * @throw std::domain_error If the input vector is size 0.
 */
template <typename Mat, require_eigen_t<Mat>* = nullptr,
          require_not_eigen_vector_t<Mat>* = nullptr,
          require_t<is_fvar<value_type_t<Mat>>>* = nullptr>
inline auto log_softmax(const Mat& m) {
  check_nonzero_size("log_softmax", "m", m);
  const auto& m_ref = to_ref(m);
  const auto val = m_ref.val().eval();
  const auto shifted
      = (val.array().colwise() - val.rowwise().maxCoeff().array()).eval();
  const auto exp_s = shifted.exp().eval();
  const auto row_sums = exp_s.rowwise().sum().eval();
  const auto lsm_val = (shifted.colwise() - row_sums.log()).matrix().eval();
  // softmax values needed for the tangent: d_in - softmax(x) ⊙ dot(softmax(x), d_in)
  const auto s = (exp_s.colwise() / row_sums).eval();
  const auto d_in = m_ref.d().eval();
  const auto dots = (s.array() * d_in.array()).rowwise().sum().eval();
  plain_type_t<Mat> result(m_ref.rows(), m_ref.cols());
  result.val() = lsm_val;
  result.d() = (d_in.array().colwise() - dots.array()).matrix();
  return result;
}

template <typename RowVec, require_eigen_row_vector_t<RowVec>* = nullptr,
          require_t<is_fvar<value_type_t<RowVec>>>* = nullptr>
inline auto log_softmax(const RowVec& x) {
  return log_softmax(x.transpose()).transpose().eval();
}

template <typename T, require_vector_st<is_fvar, T>* = nullptr,
          require_not_t<is_eigen_row_vector<std::decay_t<T>>>* = nullptr>
inline auto log_softmax(T&& x) {
  return apply_vector_unary<T>::apply(std::forward<T>(x), [](auto&& alpha) {
    using T_alpha = std::decay_t<decltype(alpha)>;
    using T_fvar = value_type_t<T_alpha>;
    using T_inner = typename T_fvar::Scalar;

    auto&& alpha_ref = to_ref(std::forward<decltype(alpha)>(alpha));
    const Eigen::Matrix<T_inner, -1, 1> val = alpha_ref.val();
    const Eigen::Matrix<T_inner, -1, 1> s = softmax(val);
    const auto d_in = alpha_ref.d().eval();
    const T_inner dot_sd = s.dot(d_in);

    Eigen::Matrix<T_fvar, -1, 1> result(alpha_ref.size());
    result.val() = log_softmax(val);
    result.d() = (d_in.array() - dot_sd).matrix();
    return result;
  });
}

}  // namespace math
}  // namespace stan
#endif
