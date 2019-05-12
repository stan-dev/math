#ifndef STAN_MATH_PRIM_MAT_META_MATRIX_RETURN_TYPE_HPP
#define STAN_MATH_PRIM_MAT_META_MATRIX_RETURN_TYPE_HPP

#include <Eigen/Core>
#include <stan/math/prim/mat/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>

namespace stan {
namespace math {

/**
 * Metaprogram that determines if the return type
 * is an Eigen::Matrix of partials_return_type or
 * the partials_return_type.
 * @tparam T Type to test
 */
template <typename... T>
using matrix_return_type = typename std::conditional<is_vector<T...>::value,
      Eigen::Matrix<typename partials_return_type<T...>::type, -1, 1>,
      typename partials_return_type<T...>::type>::type;

}  // namespace math
}  // namespace stan

#endif
