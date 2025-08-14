#ifndef STAN_MATH_PRIM_FUN_DIAGONAL_HPP
#define STAN_MATH_PRIM_FUN_DIAGONAL_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return a column vector of the diagonal elements of the
 * specified matrix.  The matrix is not required to be square.
 *
 * @tparam T type of the matrix
 * @param m Specified matrix.
 * @return Diagonal of the matrix.
 */
template <typename T, require_matrix_t<T>* = nullptr>
inline auto diagonal(T&& m) {
  return make_holder([](auto&& m_) {
    return m_.diagonal();
  }, std::forward<T>(m));
}

}  // namespace math
}  // namespace stan

#endif
