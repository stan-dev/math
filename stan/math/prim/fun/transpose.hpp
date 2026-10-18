#ifndef STAN_MATH_PRIM_FUN_TRANSPOSE_HPP
#define STAN_MATH_PRIM_FUN_TRANSPOSE_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Transposes a matrix.
 * @tparam T type of the matrix or expression
 * @param m matrix or expression
 * @return transposed matrix
 */
template <typename T, require_matrix_t<T>* = nullptr>
inline auto transpose(T&& m) {
  return make_holder([](auto&& m_) { return m_.transpose(); },
                     std::forward<T>(m));
}

}  // namespace math
}  // namespace stan
#endif
