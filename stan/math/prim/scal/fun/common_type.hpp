#ifndef STAN_MATH_PRIM_SCAL_FUN_COMMON_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_FUN_COMMON_TYPE_HPP

#include <boost/math/tools/promotion.hpp>

namespace stan {
  namespace math {
    /**
     * Struct carries out type promotion.
     *
     * @tparam T1 type of arg1
     * @tparam T2 type of arg2
     */
    template <typename T1, typename T2>
    struct common_type {
      typedef typename boost::math::tools::promote_args<T1, T2>::type type;
    };

  }
}

#endif
