#ifndef STAN_MATH_PRIM_SCAL_FUN_SELFSTART_ASYMP_HPP
#define STAN_MATH_PRIM_SCAL_FUN_SELFSTART_ASYMP_HPP

#include <vector>
#include <cmath>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>

namespace stan {

  namespace math {

    template <typename T1, typename T2, typename T3, typename T4>
    inline typename boost::math::tools::promote_args<T1, T2, T3, T4>::type
    selfstart_asymp(const T1& x, const T2& Asym, const T3& R0,
                    const T4& exp_lrc) {
      using std::exp;
      return Asym + (R0 - Asym) * exp((-exp_lrc * x));
    }	

    template <typename T1, typename T2, typename T3, typename T4> inline
    std::vector<typename boost::math::tools::promote_args<T1, T2, T3, T4>::type>
    selfstart_asymp(const std::vector<T1>& x, const T2& Asym, const T3& R0,
                    const T4& lrc) {

      using std::exp;
      using stan::math::check_finite;
      using stan::math::check_not_nan;

      static const char* function("stan::math::selfstart_asymp");
      std::vector<typename
      boost::math::tools::promote_args<T1, T2, T3, T4>::type>
      value(x.size());

      if (x.size() == 0) return value;

      check_not_nan(function, "Random variable", x);
      check_finite(function, "Asymptote variable", Asym);
      check_finite(function, "Zero Input variable", R0);
      check_finite(function, "Natural log of rate constant", lrc);

      T4 exp_lrc = exp(lrc);

      for (std::size_t i = 0; i < x.size(); ++i)
        value[i] = selfstart_asymp(x[i], Asym, R0, exp_lrc);

      return value;
    }

    template <typename T1, typename T2, typename T3, typename T4> inline
    std::vector<typename boost::math::tools::promote_args<T1, T2, T3, T4>::type>
    selfstart_offset_asymp(const std::vector<T1>& x, const T2& Asym,
                           const T3& lrc, const T4& c0) {
      using std::exp;
      using stan::math::check_finite;
      using stan::math::check_not_nan;

      static const char* function("stan::math::selfstart_offset_asymp");
      std::vector<typename
      boost::math::tools::promote_args<T1, T2, T3, T4>::type>
      value(x.size());

      if (x.size() == 0) return value;

      check_not_nan(function, "Random variable", x);
      check_finite(function, "Asymptote variable", Asym);
      check_finite(function, "Natural log of rate constant", lrc);
      check_finite(function, "Zero Response variable", c0);

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
