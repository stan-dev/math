#ifndef STAN_MATH_REV_SCAL_META_PARTIALS_TYPE_HPP
#define STAN_MATH_REV_SCAL_META_PARTIALS_TYPE_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/meta/partials_type.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <type_traits>

namespace stan {

/**
 * Specialization of partials type returns double if input type is a double.
 */
template <typename T>
struct partials_type<T, require_var_t<T>> {
  using type = typename std::decay_t<T>::Scalar;
};

}  // namespace stan
#endif
