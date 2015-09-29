#ifndef STAN_MATH_PRIM_SCAL_FUN_SELFSTART_OFFSET_ASYMP_HPP
#define STAN_MATH_PRIM_SCAL_FUN_SELFSTART_OFFSET_ASYMP_HPP

#include <vector>
#include <cmath>
#include <boost/math/tools/promotion.hpp>

namespace stan {

  namespace math {

    template <typename T1, typename T2, typename T3, typename T4> inline
    std::vector<typename boost::math::tools::promote_args<T1, T2, T3, T4>::type>
    selfstart_offset_asymp(const std::vector<T1>& x, const T2& Asym,
                           const T3& lrc, const T4& c0) {
      using std::exp;

      std::vector<typename
      boost::math::tools::promote_args<T1, T2, T3, T4>::type>
      value(x.size());

      T3 exp_lrc = exp(lrc);

      for (std::size_t i = 0; i < x.size(); ++i)
        value[i] = Asym * (1 - exp(-exp_lrc * (x[i] - c0)));

      return value;
    }

    template <typename T1, typename T2, typename T3> inline
    std::vector<typename boost::math::tools::promote_args<T1, T2, T3>::type>
    selfstart_origin_asymp(const std::vector<T1>& x, const T2& Asym,
                           const T3& lrc) {
      return selfstart_offset_asymp(x, Asym, lrc, 0);
    }


  }

}

#endif
