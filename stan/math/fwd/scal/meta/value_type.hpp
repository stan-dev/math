#ifndef STAN_MATH_FWD_SCAL_META_VALUE_TYPE_HPP
#define STAN_MATH_FWD_SCAL_META_VALUE_TYPE_HPP

#include <stan/math/fwd/scal/meta/is_fvar.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Template metaprogram class to compute the type of values stored
 * in a standard vector.
 *
 * @tparam T type of elements in standard vector.
 */
template <typename T>
struct value_type<T, std::enable_if_t<is_fvar<T>::value>> {
  using type = typename std::decay_t<T>::Scalar;
};

}  // namespace math
}  // namespace stan
#endif
