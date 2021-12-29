#ifndef STAN_MATH_PRIM_FUN_NORM1_HPP
#define STAN_MATH_PRIM_FUN_NORM1_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cstddef>
#include <vector>

namespace stan {
namespace math {

inline double norm1(const std::vector<double>& x) {
  double norm = 0.0;
  for (double i : x) {
    norm += std::abs(i);
  }
  return norm;
}

/**
 * Returns L1 norm of a vector. For vectors that equals the
 * sum of magnitudes of its individual elements.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase)
 * @param v Vector.
 * @return L1 norm of v.
 */
template <typename T, require_eigen_t<T>* = nullptr,
          require_not_eigen_vt<is_var, T>* = nullptr>
inline auto norm1(const T& v) {
  return v.template lpNorm<1>();
}

}  // namespace math
}  // namespace stan

#endif
