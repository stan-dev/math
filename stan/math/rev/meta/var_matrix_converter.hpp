#ifndef STAN_MATH_REV_META_VAR_MATRIX_CONVERTER_HPP
#define STAN_MATH_REV_META_VAR_MATRIX_CONVERTER_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/meta/is_var.hpp>
#include <stan/math/rev/core/var.hpp>

namespace stan {

namespace internal {
  template <typename T, typename = void>
  struct var_matrix_converter;
  template <typename T>
  struct var_matrix_converter<T, require_eigen_t<T>> {
    using type = math::var_value<math::promote_scalar_t<double, std::decay_t<T>>>;
  };

  template <typename T>
  struct var_matrix_converter<T, require_var_t<T>> {
    using type = T;
  };
}

template <typename T>
using var_matrix_converter_t = typename internal::var_matrix_converter<std::decay_t<T>>::type;

}

#endif
