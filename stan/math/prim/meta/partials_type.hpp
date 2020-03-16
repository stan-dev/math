#ifndef STAN_MATH_PRIM_META_PARTIALS_TYPE_HPP
#define STAN_MATH_PRIM_META_PARTIALS_TYPE_HPP

namespace stan {

/** \ingroup type_trait
 * This base implementation will contain a static member function named type
 * equal to the type passed into it. When this is specialized for vars the type
 * will be double and fvar<T> will have a member type of value T.
 */
template <typename T, typename = void>
struct partials_type {
  using type = T;
};

/** \ingroup type_trait
 * Helper alias for accessing the partial type.
 */
template <typename T>
using partials_type_t = typename partials_type<T>::type;

}  // namespace stan
#endif
