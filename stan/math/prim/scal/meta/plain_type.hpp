#ifndef STAN_MATH_PRIM_SCAL_META_PLAIN_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_META_PLAIN_TYPE_HPP

#include <type_traits>

namespace stan {
namespace math {

/**
 * Determines plain (non expression) type associated with \c T. For non \c Eigen
 * types it is the input type.
 * @tparam T type to determine plain type of
 */
template <typename T, typename Enable = void>
struct plain_type {
  using type = std::decay_t<T>;
};

template <typename T>
using plain_type_t = typename plain_type<T>::type;

}  // namespace math
}  // namespace stan

#endif  // STAN_MATH_PRIM_SCAL_META_PLAIN_TYPE_HPP
