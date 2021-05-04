#ifndef STAN_MATH_PRIM_FUN_INITIALIZE_FILL_HPP
#define STAN_MATH_PRIM_FUN_INITIALIZE_FILL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Fill the specified container with the specified value.
 *
 * The specified matrix is filled by element.
 *
 * @tparam EigMat Type inheriting from `EigenBase`
 * @tparam S Type of value.
 *
 * @param x Container.
 * @param y Value.
 */
template <typename EigMat, typename S, require_eigen_t<EigMat>* = nullptr,
          require_stan_scalar_t<S>* = nullptr>
inline void initialize_fill(EigMat& x, const S& y) {
  x.fill(y);
}

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
template <
    typename T, typename S,
    require_t<std::is_assignable<std::decay_t<T>&, std::decay_t<S>>>* = nullptr>
inline void initialize_fill(T& x, S&& y) {
  x = std::forward<S>(y);
}

/**
 * Fill the specified container with the specified value.
 *
 * Each container in the specified standard vector is filled
 * recursively by calling <code>fill</code>.
 *
 * @tparam Vec A standard vector
 * @tparam S type of value
 * @param[in] x Container.
 * @param[in, out] y Value.
 */
template <typename Vec, typename S, require_std_vector_t<Vec>* = nullptr>
inline void initialize_fill(Vec& x, S&& y) {
  for (auto& x_val : x) {
    initialize_fill(x_val, y);
  }
}

}  // namespace math
}  // namespace stan

#endif
