#ifndef STAN_MATH_PRIM_SCAL_FUN_SELFSTART_WEIBULL_HPP
#define STAN_MATH_PRIM_SCAL_FUN_SELFSTART_WEIBULL_HPP

#include <vector>
#include <cmath>
#include <boost/math/tools/promotion.hpp>

namespace stan {

  namespace math {

    template <typename T1, typename T2, typename T3, typename T4, typename T5>
    inline std::vector<typename
    boost::math::tools::promote_args<T1, T2, T3, T4, T5>::type>
    selfstart_weibull(const std::vector<T1>& x, const T2& Asym, const T3& Drop,
                      const T4& lrc, const T5& pwr) {
      using std::exp;
      using std::pow;

      std::vector<typename
      boost::math::tools::promote_args<T1, T2, T3, T4, T5>::type>
      value(x.size());

      T4 exp_lrc = exp(lrc);

      for (std::size_t i = 0; i < x.size(); ++i)
        value[i] = Asym - Drop * exp(-exp_lrc * pow(x[i], pwr));

      return value;
    }

  }

}
#endif
