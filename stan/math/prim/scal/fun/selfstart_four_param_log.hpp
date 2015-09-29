#ifndef STAN_MATH_PRIM_SCAL_FUN_SELFSTART_FOUR_PARAM_LOG_HPP
#define STAN_MATH_PRIM_SCAL_FUN_SELFSTART_FOUR_PARAM_LOG_HPP

#include <vector>
#include <cmath>
#include <boost/math/tools/promotion.hpp>

namespace stan {

  namespace math {

    template <typename T1, typename T2, typename T3, typename T4, typename T5>
    inline std::vector<typename
    boost::math::tools::promote_args<T1, T2, T3, T4, T5>::type>
    selfstart_four_param_log(const std::vector<T1>& x, const T2& A,
                             const T3& B, const T4& xmid, const T5& scal) {
      using std::exp;

      std::vector<typename
      boost::math::tools::promote_args<T1, T2, T3, T4, T5>::type>
      value(x.size());

      for (std::size_t i = 0; i < x.size(); ++i)
        value[i] = A + (B - A)/(1 + exp((xmid - x[i])/scal));

      return value;
    }
  }

}
#endif
