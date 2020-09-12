#ifndef STAN_MATH_OPENCL_IS_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_IS_MATRIX_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {

namespace math {
/**
 * Dummy class to instantiate matrix_cl to enable for specific types.
 * @ingroup matrix_cl_group
 */
template <typename T, typename = void>
class matrix_cl {
 public:
  using Scalar = T;
  using type = T;
};
}  // namespace math

namespace internal {

/** \ingroup type_traits
 * @internal
 * This underlying implementation is used when the type is not an std vector.
 */
template <typename T>
struct is_matrix_cl_impl : std::false_type {};

/** \ingroup type_traits
 * @internal
 * This specialization implementation has a static member named value when the
 * template type is an std vector.
 */
template <typename... Args>
struct is_matrix_cl_impl<stan::math::matrix_cl<Args...>> : std::true_type {};

}  // namespace internal

template <typename T, typename = void>
struct is_matrix_cl : std::false_type {};

/** \ingroup type_traits
 * Checks if the decayed type of T is a matrix_cl.
 */
template <typename T>
struct is_matrix_cl<
    T, std::enable_if_t<internal::is_matrix_cl_impl<std::decay_t<T>>::value>>
    : std::true_type {};

STAN_ADD_REQUIRE_UNARY(matrix_cl, is_matrix_cl, matrix_cl_group);
STAN_ADD_REQUIRE_CONTAINER(matrix_cl, is_matrix_cl, matrix_cl_group);

}  // namespace stan
#endif
#endif
