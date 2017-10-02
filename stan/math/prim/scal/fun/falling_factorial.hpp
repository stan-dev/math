#ifndef STAN_MATH_PRIM_SCAL_FUN_FALLING_FACTORIAL_HPP
#define STAN_MATH_PRIM_SCAL_FUN_FALLING_FACTORIAL_HPP

#include <boost/math/special_functions/factorials.hpp>
#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <limits>
#include <string>

namespace stan {
  namespace math {

    /**
     *
       \f[
       \mbox{falling\_factorial}(x, n) =
       \begin{cases}
         \textrm{error} & \mbox{if } x \leq 0\\
         (x)_n & \mbox{if } x > 0 \textrm{ and } -\infty \leq n \leq \infty \\[6pt]
         \textrm{NaN} & \mbox{if } x = \textrm{NaN or } n = \textrm{NaN}
       \end{cases}
       \f]

       \f[
       \frac{\partial\, \mbox{falling\_factorial}(x, n)}{\partial x} =
       \begin{cases}
         \textrm{error} & \mbox{if } x \leq 0\\
         \frac{\partial\, (x)_n}{\partial x} & \mbox{if } x > 0 \textrm{ and } -\infty \leq n \leq \infty \\[6pt]
         \textrm{NaN} & \mbox{if } x = \textrm{NaN or } n = \textrm{NaN}
       \end{cases}
       \f]

       \f[
       \frac{\partial\, \mbox{falling\_factorial}(x, n)}{\partial n} =
       \begin{cases}
         \textrm{error} & \mbox{if } x \leq 0\\
         \frac{\partial\, (x)_n}{\partial n} & \mbox{if } x > 0 \textrm{ and } -\infty \leq n \leq \infty \\[6pt]
         \textrm{NaN} & \mbox{if } x = \textrm{NaN or } n = \textrm{NaN}
       \end{cases}
       \f]

       \f[
       (x)_n=\frac{\Gamma(x+1)}{\Gamma(x-n+1)}
       \f]

       \f[
       \frac{\partial \, (x)_n}{\partial x} = (x)_n\Psi(x+1)
       \f]

       \f[
       \frac{\partial \, (x)_n}{\partial n} = -(x)_n\Psi(n+1)
       \f]
     *
     */
    template<typename T1, typename T2>
    inline typename boost::math::tools::promote_args<T1, T2>::type
    falling_factorial(const T1& x, const T2&n) {
      if (is_nan(x) || is_nan(n))
        return std::numeric_limits<double>::quiet_NaN();
      static const std::string function = "falling_factorial";
      check_nonnegative(function, "second argument", n);
      return boost::math::falling_factorial(x, n, boost_policy_t());
    }

  }
}

#endif
