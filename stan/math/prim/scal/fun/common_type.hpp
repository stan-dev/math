#ifndef STAN_MATH_PRIM_SCAL_FUN_COMMON_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_FUN_COMMON_TYPE_HPP

#include <boost/math/tools/promotion.hpp>

namespace stan {
  namespace math {
    /**
     * Struct which calculates wider type given two types.
     * See <a
     * href="http://boost.org/doc/libs/1_63_0/boost/math/tools/promotion.hpp">
     * Boost/math/tools/promotion.hpp</a> for definition of "wider".
     *
     * <p>This is the base implementation for scalar types.
     *
     */
    template <typename T1, typename T2>
    struct common_type {
      typedef typename boost::math::tools::promote_args<T1, T2>::type type;
    };

  }
}

#endif
