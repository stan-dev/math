#ifndef STAN_MATH_PRIM_MAT_FUN_HEAD_HPP
#define STAN_MATH_PRIM_MAT_FUN_HEAD_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/check_column_index.hpp>
#include <stan/math/prim/mat/err/check_row_index.hpp>
#include <stan/math/prim/mat/err/check_std_vector_index.hpp>
#include <stan/math/prim/mat/vectorize/apply_vector_unary.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the specified number of elements as a vector
 * from the front of the specified vector.
 *
 * @tparam T Type of input vector.
 * @tparam T2 Type of size variable.
 * @param x Vector input.
 * @param n Size of return.
 * @return The first n elements of v.
 * @throw std::out_of_range if n is out of range.
 */
template <typename T, typename T2>
inline auto head(T&& x, const T2& n) {
  return apply_vector_unary<T>::apply_scalar(
      std::forward<T>(x), n, [&](auto& v, auto& m) {
        if (m != 0) {
          if (v.rows() == 1) {
            check_column_index("head", "n", v, m);
          } else {
            check_row_index("head", "n", v, m);
          }
        }
        return v.head(m);
      });
}

}  // namespace math
}  // namespace stan
#endif
