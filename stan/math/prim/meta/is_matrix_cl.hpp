#ifndef STAN_MATH_PRIM_META_IS_MATRIX_CL_HPP
#define STAN_MATH_PRIM_META_IS_MATRIX_CL_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Non-templated base class for `matrix_cl` simplifies checking if something is
 * matrix_cl.
 */
class matrix_cl_base {};

}  // namespace math

/** \ingroup type_traits
 * Checks if the decayed type of T is a matrix_cl.
 */
template <typename T>
struct is_matrix_cl
    : public std::is_base_of<math::matrix_cl_base, std::decay_t<T>> {};

STAN_ADD_REQUIRE_UNARY(matrix_cl, is_matrix_cl, matrix_cl_group);
STAN_ADD_REQUIRE_CONTAINER(matrix_cl, is_matrix_cl, matrix_cl_group);

}  // namespace stan
#endif
