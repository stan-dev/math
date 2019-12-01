#ifndef STAN_MATH_PRIM_META_AD_PROMOTABLE_HPP
#define STAN_MATH_PRIM_META_AD_PROMOTABLE_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <type_traits>

namespace stan {
namespace math {

/** \ingroup type_trait
 * If the type From can be converted to To using implicit conversions, or
 * both From and To are possibly cv-qualified void),
 * provides the member constant value equal to true.
 *
 * @tparam From promoted type
 * @tparam To target type
 */
template <typename From, typename To, typename = void>
struct ad_promotable
    : bool_constant<
          std::is_convertible<std::decay_t<From>, std::decay_t<To>>::value> {};

}  // namespace math
}  // namespace stan
#endif
