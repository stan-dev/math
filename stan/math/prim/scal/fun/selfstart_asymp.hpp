#ifndef STAN_MATH_PRIM_SCAL_FUN_SELFSTART_ASYMP_HPP
#define STAN_MATH_PRIM_SCAL_FUN_SELFSTART_ASYMP_HPP

#include <vector>
#include <boost/math/tools/promotion.hpp>

namespace stan {
  
  namespace math {

    template <typename T1, typename T2, typename T3> 
    inline typename boost::math::tools::promote_args<T1, T2, T3>::type 
    selfstart_asymp(double x, const T1& Asym, const T2& R0, const T3& exp_lrc) {
      using std::exp;
      return Asym + (R0 - Asym) * exp((-exp_lrc * x));
    }	
    
    template <typename T1, typename T2, typename T3> inline
    std::vector<typename boost::math::tools::promote_args<T1, T2, T3>::type> 
    selfstart_asymp(const std::vector<double>& x, const T1& Asym,
                    const T2& R0, const T3& lrc) {
      using std::exp;

      std::vector<typename boost::math::tools::promote_args<T1, T2, T3>::type> 
      value(x.size());
      
      T3 exp_lrc = exp(lrc);

      for (std::size_t i = 0; i < x.size(); ++i)
        value[i] = selfstart_asymp(x[i], Asym, R0, exp_lrc);

      return value; 
    }

  }
}

#endif
