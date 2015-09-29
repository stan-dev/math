#ifndef STAN_MATH_PRIM_SCAL_FUN_SELFSTART_FIRST_ORDER_COMPART_HPP
#define STAN_MATH_PRIM_SCAL_FUN_SELFSTART_FIRST_ORDER_COMPART_HPP

#include <vector>
#include <cmath>
#include <boost/math/tools/promotion.hpp>

namespace stan {

  namespace math {

    template <typename T1, typename T2, typename T3, typename T4, typename T5>
    inline std::vector<typename
    boost::math::tools::promote_args<T1, T2, T3, T4, T5>::type>
    selfstart_first_order_compart(const std::vector<T1>& x, const T2& Dose,
                                  const T3& lKe, const T4& lKa,
                                  const T5& lCl) {

      using std::exp;

      std::vector<typename
      boost::math::tools::promote_args<T1, T2, T3, T4, T5>::type>
      value(x.size());

      T3 exp_lke = exp(lKe);
      T4 exp_lka = exp(lKa);
      typename boost::math::tools::promote_args<T2, T3, T4, T5>::type
      dose_const = Dose * exp_lke * exp_lka/(exp(lCl) * (exp_lka - exp_lke));

      for (std::size_t i = 0; i < x.size(); ++i)
        value[i] = dose_const * (exp(-exp_lke * x[i]) - exp(-exp_lka * x[i]));

      return value;
    }
  }

}
#endif
