#ifndef STAN_MATH_PRIM_FUN_DOT_SELF_HPP
#define STAN_MATH_PRIM_FUN_DOT_SELF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cstddef>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the dot product of the specified vector with itself.
 *
 * @tparam StdVec a standard vector.
 * @param v Vector.
 * @throw std::domain_error If v is not vector dimensioned.
 */
template <typename StdVec, require_std_vector_t<StdVec>* = nullptr>
inline auto dot_self(StdVec&& x) {
  value_type_t<StdVec> sum = 0.0;
  for (auto&& i : x) {
    sum += i * i;
  }
  return sum;
}

/**
 * Returns the dot product of the specified vector with itself.
 *
 * @tparam EigVec A type deriving from `Eigen::MatrixBase` with compile time
 *  rows or columns equal to 1.
 * @param v Vector.
 * @throw std::domain_error If v is not vector dimensioned.
 */
template <typename EigMat,
          require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr>
inline double dot_self(EigMat&& v) {
  return v.squaredNorm();
}

}  // namespace math
}  // namespace stan

#endif
