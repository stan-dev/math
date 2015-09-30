#ifndef STAN_MATH_PRIM_SCAL_FUN_SELFSTART_MICMEN_HPP
#define STAN_MATH_PRIM_SCAL_FUN_SELFSTART_MICMEN_HPP

#include <vector>
#include <cmath>
#include <boost/math/tools/promotion.hpp>

namespace stan {

  namespace math {

    template <typename T1, typename T2, typename T3>
    inline std::vector<typename
    boost::math::tools::promote_args<T1, T2, T3>::type>
    selfstart_micmen(const std::vector<T1>& x, const T2& Vm, const T3& K) {
      using std::exp;

      std::vector<typename boost::math::tools::promote_args<T1, T2, T3>::type>
      value(x.size());

      for (std::size_t i = 0; i < x.size(); ++i)
        value[i] = Vm * x[i]/(K + x[i]);

      return value;
    }

  }

}
#endif
