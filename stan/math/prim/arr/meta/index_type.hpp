#ifndef STAN_MATH_PRIM_ARR_META_INDEX_TYPE_HPP
#define STAN_MATH_PRIM_ARR_META_INDEX_TYPE_HPP

#include <stan/math/prim/arr/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/index_type.hpp>
#include <stan/math/prim/scal/meta/require_generics.hpp>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

/**
 * Template metaprogram class to compute the type of index for a
 * standard vector.
 *
 * @tparam T type of elements in standard vector.
 */
template <typename T>
struct index_type<T, require_std_vector_t<T>> {
  using type = typename std::decay_t<T>::size_type;
};

}  // namespace math
}  // namespace stan

#endif
