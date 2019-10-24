#ifndef STAN_MATH_PRIM_ARR_META_IS_CONSTANT_HPP
#define STAN_MATH_PRIM_ARR_META_IS_CONSTANT_HPP

#include <stan/math/prim/arr/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/bool_constant.hpp>
#include <stan/math/prim/scal/meta/is_constant.hpp>
#include <stan/math/prim/scal/meta/require_generics.hpp>
#include <type_traits>
#include <vector>

namespace stan {
/**
 * Defines a static member named value and sets it to true
 * if the type of the elements in the provided std::vector
 * is constant, false otherwise. This is used in
 * the is_constant_all metaprogram.
 * @tparam type of the elements in the std::vector
 */
template <typename T>
struct is_constant<T, require_std_vector_t<T>>
    : bool_constant<is_constant<typename std::decay_t<T>::value_type>::value> {
};

}  // namespace stan
#endif
