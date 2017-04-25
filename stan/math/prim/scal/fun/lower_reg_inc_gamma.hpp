#ifndef STAN_MATH_PRIM_SCAL_FUN_LOWER_REG_INC_GAMMA_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LOWER_REG_INC_GAMMA_HPP

#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/gamma_p.hpp>
#include <stan/math/prim/scal/fun/inv.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/prim/scal/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>

namespace stan {
  namespace math {

    template <int IMPL=0>
    struct lower_reg_inc_gamma {
      template <int ARG, typename T1, typename T2>
      static 
      typename boost::math::tools::promote_args<T1, T2>::type
      grad(T1 a, T2 z, double precision = 1e-10, int max_steps = 1e5) {
        typedef typename boost::math::tools::promote_args<T1, T2>::type TP;
        TP d = std::numeric_limits<TP>::quiet_NaN;
        return d;
      }
    };
     

    template <>
    struct lower_reg_inc_gamma<1> {
      template <int ARG, typename T1, typename T2>
      static
      typename boost::math::tools::promote_args<T1, T2>::type
      grad(T1 a, T2 z, double precision = 1e-10,
                                 int max_steps = 1e5) {
        typedef typename boost::math::tools::promote_args<T1, T2>::type TP;
  
        if (is_nan(a)) return std::numeric_limits<TP>::quiet_NaN();
        if (is_nan(z)) return std::numeric_limits<TP>::quiet_NaN();
  
        using std::exp;
        using std::log;
  
        T2 log_z = log(z);
  
        int k = 0;
        TP sum_a = 0.0;
        TP sum_b = 0.0;
  
        TP a_plus_n = a;  // n begins at zero.
        TP lgamma_a_plus_n_plus_1 = lgamma(a+1);
        TP term;
        while (true) {
          term = exp(a_plus_n * log_z - lgamma_a_plus_n_plus_1); 
          sum_a += term;
          if (term <= precision)
            break;  
          if (k >= max_steps)
            domain_error("d_lower_reg_inc_gamma_da", "k (internal counter)", 
              max_steps, "exceeded ", " iterations, gamma_p(a,z) gradient (a) "
              "did not converge.");
          k++;
          lgamma_a_plus_n_plus_1 += log(a_plus_n + 1);
          a_plus_n++;
        }
  
        k = 0;
        a_plus_n = a;  // n begins at zero.
        lgamma_a_plus_n_plus_1 = lgamma(a+1);
        while (true) {
          term = exp(a_plus_n * log_z - lgamma_a_plus_n_plus_1 + 
            log(digamma(a_plus_n + 1)));
          sum_b += term;
          if (term <= precision)
            return (log_z/exp(z)) * sum_a - inv(exp(z)) * sum_b;
          if (k >= max_steps)
            domain_error("d_lower_reg_inc_gamma_da", "k (internal counter)", 
              max_steps, "exceeded ", " iterations, gamma_p(a,z) gradient (a) "
              "did not converge.");
          k++;
          lgamma_a_plus_n_plus_1 += log(a_plus_n + 1);
          a_plus_n++;
        }
      }
    };
  
  
    template <>
    struct lower_reg_inc_gamma<2> {
      template <int ARG, typename T1, typename T2>
      static
      typename boost::math::tools::promote_args<T1, T2>::type
      grad(T1 a, T2 z, double precision = 1e-10,
                                 int max_steps = 1e5) {
        typedef typename boost::math::tools::promote_args<T1, T2>::type TP;
  
        if (is_nan(a)) return std::numeric_limits<TP>::quiet_NaN();
        if (is_nan(z)) return std::numeric_limits<TP>::quiet_NaN();
        double tg = tgamma(a);
        double dig = digamma(a); 
        return -stan::math::grad_reg_inc_gamma(a, z, tg, dig, max_steps, precision);                   
      }
    };

    template <>
    struct lower_reg_inc_gamma<3> {
      template <int ARG, typename T1, typename T2>
      static
      typename boost::math::tools::promote_args<T1, T2>::type
      grad(T1 a, T2 z, double precision = 1e-10,
                                 int max_steps = 1e5) {
        typedef typename boost::math::tools::promote_args<T1, T2>::type TP;
  
        if (is_nan(a)) return std::numeric_limits<TP>::quiet_NaN();
        if (is_nan(z)) return std::numeric_limits<TP>::quiet_NaN();
  
        using std::exp;
        using std::log;
        using std::pow;
  
        T2 log_z = log(z);
  
        int k = 0;
        TP sum_a = 0.0;
  
        TP a_plus_n = a;  // n begins at zero.
        TP lgamma_a_plus_n_plus_1 = lgamma(a+1);
        TP term;
        while (true) {
          term = exp(a_plus_n * log_z - lgamma_a_plus_n_plus_1); 
          sum_a += term;
          if (term <= precision)
            break;  
          if (k >= max_steps)
            domain_error("d_lower_reg_inc_gamma_da", "k (internal counter)", 
              max_steps, "exceeded ", " iterations, gamma_p(a,z) gradient (a) "
              "did not converge.");
          k++;
          lgamma_a_plus_n_plus_1 += log(a_plus_n + 1);
          a_plus_n++;
        }
  
        int n = 1;
        TP sum_b = pow(z,a) * digamma(a + 1) / tgamma(a + 1);
        a_plus_n = a + 1;  // n begins at one.
        lgamma_a_plus_n_plus_1 = lgamma(a+2);
        while (true) {
          term = exp(a_plus_n * log_z - lgamma_a_plus_n_plus_1 + 
            log(digamma(a_plus_n + 1)));
          sum_b += term;
          if (term <= precision)
            return (log_z/exp(z)) * sum_a - inv(exp(z)) * sum_b;
          if (n >= max_steps)
            domain_error("d_lower_reg_inc_gamma_da", "n (internal counter)", 
              max_steps, "exceeded ", " iterations, gamma_p(a,z) gradient (a) "
              "did not converge.");
          n++;
          lgamma_a_plus_n_plus_1 += log(a_plus_n + 1);
          a_plus_n++;
        }
      }
    };

    template <>
    struct lower_reg_inc_gamma<0> {
      template <int ARG, typename T1, typename T2>
      static
      typename boost::math::tools::promote_args<T1, T2>::type
      grad(T1 a, T2 z, double precision = 1e-10,
                                 int max_steps = 1e5) {
        typedef typename boost::math::tools::promote_args<T1, T2>::type TP;
        TP result;
        if ((a < 0.8 && z > 15.0) || (a < 12.0 && z > 30.0) ||
            a < sqrt(-756 - z * z + 60 * z)) {
          result = stan::math::lower_reg_inc_gamma<2>::grad<1>(a, z, precision, max_steps);
        } else if (a < 0.8 && z <= 15.0) {
          result = stan::math::lower_reg_inc_gamma<3>::grad<1>(a, z, precision, max_steps);
        } else {
          result = stan::math::lower_reg_inc_gamma<1>::grad<1>(a, z, precision, max_steps);
        }
        return result;
      }
    };

  }
}
  

#endif
