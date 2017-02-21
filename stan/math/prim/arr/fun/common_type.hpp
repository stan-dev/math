#ifndef STAN_MATH_PRIM_ARR_FUN_COMMON_TYPE_HPP
#define STAN_MATH_PRIM_ARR_FUN_COMMON_TYPE_HPP

#include <stan/math/prim/scal/fun/common_type.hpp>
#include <boost/math/tools/promotion.hpp>
#include <vector>

namespace stan {
  namespace math {
    /**
     * Struct carries out type promotion.
     * This specialization is for vectors of differing types.
     *
     * @tparam T1 type of arg1
     * @tparam T2 type of arg2
     */
    template <typename T1, typename T2>
    struct common_type<std::vector<T1>, std::vector<T2> > {
      typedef std::vector<typename common_type<T1, T2>::type> type;
    };

  }
}

#endif
