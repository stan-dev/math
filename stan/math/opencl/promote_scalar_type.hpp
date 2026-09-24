#ifndef STAN_MATH_OPENCL_PROMOTE_SCALAR_TYPE_HPP
#define STAN_MATH_OPENCL_PROMOTE_SCALAR_TYPE_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <type_traits>

namespace stan {
namespace math {

/** \ingroup type_traits
 * Promotes the scalar type of a `matrix_cl` to an arithmetic type.
 *
 * @tparam T arithmetic type to promote to
 * @tparam S `matrix_cl` type
 */
template <typename T, typename S>
struct promote_scalar_type<
    T, S, require_all_t<std::is_arithmetic<T>, is_matrix_cl<S>>> {
  using type = matrix_cl<T>;
};

}  // namespace math
}  // namespace stan

#endif
#endif
