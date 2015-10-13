#ifndef STAN_MATH_PRIM_MAT_FUN_SELFSTART_ASYMP_HPP
#define STAN_MATH_PRIM_MAT_FUN_SELFSTART_ASYMP_HPP

#include <vector>
#include <cmath>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>

namespace stan {

  namespace math {

    /**
     * Returns the asymptotic regression function values for a vector of data
     * points x and the parameters Asym, R0, and lrc.
     *
     * This function calculates values using the following formula.
     *
     * \f$ y = Asym + (R0 - Asym) * \exp(-\exp(lrc) * x))\f$
     *
     * @param Asym horizontal asymptote as \f$ x \to \infty\f$.
     * @param R0 y intercept.
     * @param lrc natural logarithm of the rate constant.
     * @return a vector of the asymptotic regression function values
     *  for each entry in x.
     * @throw std::domain_error if any parameters are not finite or NaN.
     */
    template <typename T1, typename T2, typename T3, typename T4> inline
    std::vector<typename boost::math::tools::promote_args<T1, T2, T3, T4>::type>
    selfstart_asymp(const std::vector<T1>& x, const T2& Asym, const T3& R0,
                    const T4& lrc) {

      using std::exp;

      static const char* function("stan::math::selfstart_asymp");
      std::vector<typename
              boost::math::tools::promote_args<T1, T2, T3, T4>::type>
                       value(x.size());

      if (x.size() == 0) return value;
      check_not_nan(function, "Random variable", x);
      check_finite(function, "Random variable", x);
      check_not_nan(function, "Asymptote variable", Asym);
      check_finite(function, "Asymptote variable", Asym);
      check_not_nan(function, "Zero Input variable", R0);
      check_finite(function, "Zero Input variable", R0);
      check_not_nan(function, "Natural log of rate constant", lrc);
      check_finite(function, "Natural log of rate constant", lrc);

      T4 exp_lrc = exp(lrc);

      for (std::size_t i = 0; i < x.size(); ++i)
        value[i] = Asym + (R0 - Asym) * exp((-exp_lrc * x[i]));

      return value;
    }

    template <typename T1, typename T2, typename T3, typename T4> inline
    Eigen::Matrix<typename
             boost::math::tools::promote_args<T1, T2, T3, T4>::type,
                      Eigen::Dynamic, 1>
    selfstart_asymp(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
                       const T2& Asym, const T3& R0, const T4& lrc) {

      using std::exp;

      static const char* function("stan::math::selfstart_asymp");
      if (x.size() == 0) 
        return Eigen::Matrix<typename
                        boost::math::tools::promote_args<T1, T2, T3, T4>::type,
                                 Eigen::Dynamic, 1>();
      check_not_nan(function, "Random variable", x);
      check_finite(function, "Random variable", x);
      check_not_nan(function, "Asymptote variable", Asym);
      check_finite(function, "Asymptote variable", Asym);
      check_not_nan(function, "Zero Input variable", R0);
      check_finite(function, "Zero Input variable", R0);
      check_not_nan(function, "Natural log of rate constant", lrc);
      check_finite(function, "Natural log of rate constant", lrc);

      T4 exp_lrc = exp(lrc);

      return Asym + (R0 - Asym) * exp(-exp_lrc * x);
    }
  }

}
#endif
