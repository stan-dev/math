#ifndef STAN_MATH_FWD_FUN_SOFTMAX_HPP
#define STAN_MATH_FWD_FUN_SOFTMAX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/value_of.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/softmax.hpp>

namespace stan {
namespace math {

template <typename Mat, require_eigen_t<Mat>* = nullptr,
          require_not_eigen_vector_t<Mat>* = nullptr,
          require_t<is_fvar<value_type_t<Mat>>>* = nullptr>
inline auto softmax(const Mat& m) {
  const auto& m_ref = to_ref(m);
  const auto s = softmax(m_ref.val());
  const auto d_in = m_ref.d().eval();
  // d/dx softmax(x) applied to tangent: s ⊙ (d_in - s · d_in)  (per row)
  const auto dots = (s.array() * d_in.array()).rowwise().sum().eval();
  plain_type_t<Mat> result(m_ref.rows(), m_ref.cols());
  result.val() = s;
  result.d() = (s.array() * (d_in.array().colwise() - dots.array())).matrix();
  return result;
}

template <typename RowVec, require_eigen_row_vector_t<RowVec>* = nullptr,
          require_t<is_fvar<value_type_t<RowVec>>>* = nullptr>
inline auto softmax(const RowVec& alpha) {
  return softmax(alpha.transpose()).transpose().eval();
}

template <typename ColVec,
          require_eigen_col_vector_vt<is_fvar, ColVec>* = nullptr>
inline auto softmax(const ColVec& alpha) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using T = typename value_type_t<ColVec>::Scalar;
  if (alpha.size() == 0) {
    return Matrix<fvar<T>, Dynamic, 1>();
  }
  const auto& alpha_ref = to_ref(alpha);
  const Matrix<T, Dynamic, 1> s = softmax(value_of(alpha_ref));
  const auto d_in = alpha_ref.d().eval();
  const T dot_sd = s.dot(d_in);
  Matrix<fvar<T>, Dynamic, 1> result(alpha.size());
  result.val() = s;
  result.d() = (s.array() * (d_in.array() - dot_sd)).matrix();
  return result;
}

}  // namespace math
}  // namespace stan
#endif
