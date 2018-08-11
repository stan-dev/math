#ifndef STAN_MATH_PRIM_SCAL_META_PARTIALS_RETURN_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_META_PARTIALS_RETURN_TYPE_HPP

#include <stan/math/prim/scal/meta/partials_type.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <boost/math/tools/promotion.hpp>

namespace stan {

template <typename T, typename... T_pack>
struct partials_return_type {
  typedef typename boost::math::tools::promote_args<double,
      typename partials_type<typename scalar_type<T>::type>::type,
      typename partials_return_type<T_pack...>::type >::type type;
};

template <typename T>
struct partials_return_type<T> {
  typedef typename boost::math::tools::promote_args<double,
      typename partials_type<typename scalar_type<T>::type>::type >::type type;
};

}  // namespace stan
#endif
