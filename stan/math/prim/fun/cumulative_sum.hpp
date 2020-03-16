#ifndef STAN_MATH_PRIM_FUN_CUMULATIVE_SUM_HPP
#define STAN_MATH_PRIM_FUN_CUMULATIVE_SUM_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>
#include <numeric>
#include <functional>

namespace stan {
namespace math {

/**
 * Return the cumulative sum of the specified vector.
 *
 * The cumulative sum of a vector of values \code{x} is the
 *
 * \code x[0], x[1] + x[2], ..., x[1] + , ..., + x[x.size()-1] @endcode
 *
 * @tparam T type of elements in the vector
 * @param x Vector of values.
 * @return Cumulative sum of values.
 */
template <typename T>
inline std::vector<T> cumulative_sum(const std::vector<T>& x) {
  std::vector<T> result(x.size());
  if (x.size() == 0) {
    return result;
  }
  std::partial_sum(x.begin(), x.end(), result.begin(), std::plus<T>());
  return result;
}

/**
 * Return the cumulative sum of the specified vector.
 *
 * The cumulative sum is of the same type as the input and
 * has values defined by
 *
 * \code x(0), x(1) + x(2), ..., x(1) + , ..., + x(x.size()-1) @endcode
 *
 * @tparam EigVec type of the vector (must be derived from \c Eigen::MatrixBase
 * and have one compile time dimension equal to 1)
 *
 * @param m Vector of values.
 * @return Cumulative sum of values.
 */
template <typename EigVec, require_eigen_vector_t<EigVec>* = nullptr>
inline auto cumulative_sum(const EigVec& m) {
  using T_scalar = value_type_t<EigVec>;
  Eigen::Matrix<T_scalar, EigVec::RowsAtCompileTime, EigVec::ColsAtCompileTime>
      result(m.rows(), m.cols());
  if (m.size() == 0) {
    return result;
  }
  const Eigen::Ref<const plain_type_t<EigVec>>& m_ref = m;
  std::partial_sum(m_ref.data(), m_ref.data() + m_ref.size(), result.data(),
                   std::plus<T_scalar>());
  return result;
}

}  // namespace math
}  // namespace stan

#endif
