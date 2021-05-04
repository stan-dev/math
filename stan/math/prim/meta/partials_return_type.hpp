#ifndef STAN_MATH_PRIM_META_PARTIALS_RETURN_TYPE_HPP
#define STAN_MATH_PRIM_META_PARTIALS_RETURN_TYPE_HPP

#include <stan/math/prim/meta/partials_type.hpp>
#include <stan/math/prim/meta/promote_args.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <type_traits>

namespace stan {

/** \ingroup type_trait
 * Template metaprogram to calculate the partial derivative type resulting
 * from promoting all the scalar types of the template parameters. The
 * metaprogram can take an arbitrary number of template parameters.
 *
 * All C++ primitive types (except <code>long double</code>) are automatically
 * promoted to <code>double</code>.
 *
 * <code>partials_return_type<...></code> is a class defining a single public
 * typedef <code>type</code> that is <code>var</code> if there is a forward
 * mode variable type and is <code>double</code> otherwise (this is the most
 * common case).
 * Example usage:
 *
 *  - <code>return_type<int,double,var>::type</code> is <code>double</code>
 *  - The same thing with <code>var</code> replaced with a forward mode type
 *  like <code>fvar<T></code> will return <code>T</code>.
 *
 * @tparam T (required) A type
 * @tparam T_pack (optional) A parameter pack containing further types.
 */
template <typename T, typename... T_pack>
struct partials_return_type {
  using type = promote_args_t<double, partials_type_t<scalar_type_t<T>>,
                              typename partials_return_type<T_pack...>::type>;
};

template <typename T>
struct partials_return_type<T> {
  using type = promote_args_t<double, partials_type_t<scalar_type_t<T>>>;
};

template <typename... Args>
using partials_return_t = typename partials_return_type<Args...>::type;

}  // namespace stan
#endif
