#ifndef STAN_MATH_PRIM_SCAL_FUN_SELFSTART_BIEXP_HPP
#define STAN_MATH_PRIM_SCAL_FUN_SELFSTART_BIEXP_HPP

#include <vector>
#include <cmath>
#include <boost/math/tools/promotion.hpp>

namespace stan {

  namespace math {

    template <typename T1, typename T2, typename T3, typename T4, typename T5>
    inline std::vector<typename
    boost::math::tools::promote_args<T1, T2, T3, T4, T5>::type>
    selfstart_biexp(const std::vector<T1>& x, const T2& A1, const T3& lrc1,
                    const T4& A2, const T5& lrc2) {

      using std::exp;

      std::vector<typename
      boost::math::tools::promote_args<T1, T2, T3, T4, T5>::type>
      value(x.size());

      T3 exp_lrc1 = exp(lrc1);
      T5 exp_lrc2 = exp(lrc2);

      for (std::size_t i = 0; i < x.size(); ++i)
        value[i] = A1 * exp(-exp_lrc1 * x[i]) + A2 * exp(-exp_lrc2 * x[i]);

      return value;
    }
  }

}
#endif
