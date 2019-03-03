#ifndef STAN_MATH_PRIM_FUN_FILL_HPP
#define STAN_MATH_PRIM_FUN_FILL_HPP

#include <stan/math/prim/fun/fill.hpp>
#include <vector>
#include <stan/math/prim/fun/Eigen.hpp>


namespace stan {
namespace math {

/**
 * Fill the specified container with the specified value.
 *
 * This base case simply assigns the value to the container.
 *
 * @tparam T Type of reference container.
 * @tparam S Type of value.
 * @param x Container.
 * @param y Value.
 */
template <typename T, typename S>
void fill(T& x, const S& y) {
  x = y;
}

}  // namespace math
}  // namespace stan







namespace stan {
namespace math {

/**
 * Fill the specified container with the specified value.
 *
 * Each container in the specified standard vector is filled
 * recursively by calling <code>fill</code>.
 *
 * @tparam T Type of container in vector.
 * @tparam S Type of value.
 * @param[in] x Container.
 * @param[in, out] y Value.
 */
template <typename T, typename S>
void fill(std::vector<T>& x, const S& y) {
  for (typename std::vector<T>::size_type i = 0; i < x.size(); ++i)
    fill(x[i], y);
}

}  // namespace math
}  // namespace stan






namespace stan {
namespace math {

/**
 * Fill the specified container with the specified value.
 *
 * The specified matrix is filled by element.
 *
 * @tparam T Type of scalar for matrix container.
 * @tparam R Row type of matrix.
 * @tparam C Column type of matrix.
 * @tparam S Type of value.
 * @param x Container.
 * @param y Value.
 */
template <typename T, int R, int C, typename S>
void fill(Eigen::Matrix<T, R, C>& x, const S& y) {
  x.fill(y);
}

}  // namespace math
}  // namespace stan
#endif
