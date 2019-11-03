#ifndef STAN_MATH_PRIM_ARR_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_ARR_META_SCALAR_TYPE_HPP

#include <stan/math/prim/arr/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/require_generics.hpp>
#include <type_traits>
#include <vector>

namespace stan {
/**
 * Specialization of scalar_type for vector to recursivly return the inner
 * scalar type.
 */
template <typename T>
struct scalar_type<T, require_std_vector_t<T>> {
  using type = scalar_type_t<typename std::decay_t<T>::value_type>;
};

}  // namespace stan
#endif
