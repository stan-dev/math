#ifndef STAN_MATH_REV_META_CONDITIONAL_VAR_EIGEN_HPP
#define STAN_MATH_REV_META_CONDITIONAL_VAR_EIGEN_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/meta/plain_type.hpp>

namespace stan {

/**
 * Constructs a prim type or var_value from a scalar and a container.
 * @tparam T_scalar Determines the scalar (var/double) of the type.
 * @tparam T_container Determines the container (matrix/vector/matrix_cl ...) of
 * the type. This must be a prim type.
 */
template <typename T_scalar, typename T_container, typename = void>
struct conditional_var_value {
  using type = std::conditional_t<is_var<scalar_type_t<T_scalar>>::value,
                                  math::var_value<plain_type_t<T_container>>,
                                  plain_type_t<T_container>>;
};
template <typename T_scalar, typename T_container>
struct conditional_var_value<T_scalar, T_container,
                             require_std_vector_t<T_container>> {
  using type = std::vector<typename conditional_var_value<
      T_scalar, value_type_t<T_container>>::type>;
};

template <typename T_scalar, typename T_container>
using conditional_var_value_t =
    typename conditional_var_value<T_scalar, T_container>::type;

}  // namespace stan

#endif
