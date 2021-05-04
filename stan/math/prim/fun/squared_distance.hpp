#ifndef STAN_MATH_PRIM_FUN_SQUARED_DISTANCE_HPP
#define STAN_MATH_PRIM_FUN_SQUARED_DISTANCE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/square.hpp>

namespace stan {
namespace math {

/**
 * Returns the squared distance.
 *
 * @tparam Scal1 Type of the first scalar.
 * @tparam Scal2 Type of the second scalar.
 * @param x1 First scalar.
 * @param x2 Second scalar.
 * @return Squared distance between scalars
 * @throw std::domain_error Any scalar is not finite.
 */
template <typename Scal1, typename Scal2,
          require_all_stan_scalar_t<Scal1, Scal2>* = nullptr,
          require_all_not_var_t<Scal1, Scal2>* = nullptr>
inline return_type_t<Scal1, Scal2> squared_distance(const Scal1& x1,
                                                    const Scal2& x2) {
  check_finite("squared_distance", "x1", x1);
  check_finite("squared_distance", "x2", x2);
  return square(x1 - x2);
}

/**
 * Returns the squared distance between the specified vectors
 * of the same dimensions.
 *
 * @tparam EigVec1 type of the first vector (must be derived from \c
 * Eigen::MatrixBase and have one compile time dimension equal to 1)
 * @tparam EigVec2 type of the second vector (must be derived from \c
 * Eigen::MatrixBase and have one compile time dimension equal to 1)
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Square of distance between vectors.
 * @throw std::domain_error If the vectors are not the same
 * size.
 */
template <typename EigVec1, typename EigVec2,
          require_all_eigen_vector_t<EigVec1, EigVec2>* = nullptr,
          require_all_not_eigen_vt<is_var, EigVec1, EigVec2>* = nullptr>
inline return_type_t<EigVec1, EigVec2> squared_distance(const EigVec1& v1,
                                                        const EigVec2& v2) {
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  return (as_column_vector_or_scalar(v1) - as_column_vector_or_scalar(v2))
      .squaredNorm();
}

}  // namespace math
}  // namespace stan

#endif
