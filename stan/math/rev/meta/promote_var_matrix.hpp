#ifndef STAN_MATH_REV_META_PROMOTE_VAR_MATRIX
#define STAN_MATH_REV_META_PROMOTE_VAR_MATRIX

#include <stan/math/rev/meta/is_var.hpp>
#include <stan/math/prim/meta.hpp>


namespace stan {
  template <typename ReturnType, typename... Types>
  using promote_var_matrix_t = std::conditional_t<is_any_var_matrix<Types...>::value,
    stan::math::var_value<stan::math::promote_scalar_t<double, plain_type_t<ReturnType>>>,
     stan::math::promote_scalar_t<stan::math::var, plain_type_t<ReturnType>>>;
}

#endif
