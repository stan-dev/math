#ifndef STAN_MATH_PRIM_SCAL_META_RETURN_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_META_RETURN_TYPE_HPP

#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <boost/math/tools/promotion.hpp>

namespace stan {

/**
 * Metaprogram to calculate the base scalar return type resulting
 * from promoting all the scalar types of the template parameters.
 * Note: Sub-double types are promoted to double.
 */

template <typename T, typename... Types_pack>
struct return_type {
  typedef typename boost::math::tools::promote_args<
      double, typename scalar_type<T>::type,
      typename return_type<Types_pack...>::type>::type type;
};

template <typename T>
struct return_type<T> {
  typedef typename boost::math::tools::promote_args<
      double, typename scalar_type<T>::type>::type type;
};

}  // namespace stan
#endif
