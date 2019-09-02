#ifndef STAN_MATH_PRIM_ARR_META_VALUE_TYPE_HPP
#define STAN_MATH_PRIM_ARR_META_VALUE_TYPE_HPP

#include <stan/math/prim/arr/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/value_type.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Template metaprogram class to compute the type of values stored
 * in a standard vector.
 *
 * @tparam T type of elements in standard vector.
 */
template <typename T>
struct value_type<T, std::enable_if_t<is_std_vector<T>::value>> {
  using type = typename std::decay_t<T>::value_type;
};

}  // namespace math
}  // namespace stan
#endif
