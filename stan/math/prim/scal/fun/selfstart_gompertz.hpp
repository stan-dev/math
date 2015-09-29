#ifndef STAN_MATH_PRIM_SCAL_FUN_SELFSTART_GOMPERTZ_HPP
#define STAN_MATH_PRIM_SCAL_FUN_SELFSTART_GOMPERTZ_HPP

#include <vector>
#include <cmath>
#include <boost/math/tools/promotion.hpp>

namespace stan {

  namespace math {

    template <typename T1, typename T2, typename T3, typename T4>
    inline std::vector<typename
    boost::math::tools::promote_args<T1, T2, T3, T4>::type>
    selfstart_gompertz(const std::vector<T1>& x, const T2& Asym, const T3& b2, 
                       const T4& b3) {
      using std::exp;
      using std::pow;

      std::vector<typename
      boost::math::tools::promote_args<T1, T2, T3, T4>::type>
      value(x.size());

      for (std::size_t i = 0; i < x.size(); ++i)
        value[i] = Asym * exp(-b2 * pow(b3, x[i]));

      return value;
    }
  }

}
#endif
