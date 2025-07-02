#ifndef STAN_MATH_PRIM_FUN_SYMMETRIZE_FROM_UPPER_TRI_HPP
#define STAN_MATH_PRIM_FUN_SYMMETRIZE_FROM_UPPER_TRI_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a symmetric matrix using elements from the upper triangular part of
 * the input matrix.
 *
 * @tparam T type of elements in the matrix
 * @param m Matrix.
 * @throw std:invalid_argument if the matrix is not square.
 */
template <typename T, require_eigen_t<T>* = nullptr>
inline auto symmetrize_from_upper_tri(T&& m) {
  check_square("symmetrize_from_upper_tri", "m", m);
  return make_holder(
      [](auto&& m_) {
        return std::forward<decltype(m_)>(m_)
            .template selfadjointView<Eigen::Upper>();
      },
      std::forward<T>(m));
}

}  // namespace math
}  // namespace stan

#endif
