#ifndef STAN_MATH_FWD_SCAL_FUN_GAMMA_P_HPP
#define STAN_MATH_FWD_SCAL_FUN_GAMMA_P_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/scal/fun/gamma_p.hpp>
#include <stan/math/prim/scal/fun/grad_reg_lower_inc_gamma.hpp>

namespace stan {
  namespace math {

    template <typename T>
    inline
    fvar<T>
    gamma_p(const fvar<T>& x1, const fvar<T>& x2) {
      using std::log;
      using std::exp;
      using std::pow;
      using std::fabs;
      using boost::math::lgamma;
      using boost::math::digamma;

      T u = gamma_p(x1.val_, x2.val_);

      T der1 = grad_reg_lower_inc_gamma(x1.val_, x2.val_, 1.0e-10);
      T der2 = exp(-x2.val_ + (x1.val_ - 1.0) * log(x2.val_) - lgamma(x1.val_));

      return fvar<T>(u, x1.d_ * der1 + x2.d_ * der2);
    }

    template <typename T>
    inline
    fvar<T>
    gamma_p(const fvar<T>& x1, double x2) {
      using std::log;
      using std::exp;
      using std::pow;
      using std::fabs;
      using boost::math::tgamma;
      using boost::math::digamma;

      T u = gamma_p(x1.val_, x2);

      T der1 = grad_reg_lower_inc_gamma(x1.val_, x2, 1.0e-10);

      return fvar<T>(u, x1.d_ * der1);
    }

    template <typename T>
    inline
    fvar<T>
    gamma_p(double x1, const fvar<T>& x2) {
      using std::exp;
      using std::pow;

      T u = gamma_p(x1, x2.val_);

      double g = boost::math::tgamma(x1);

      T der2 = exp(-x2.val_) * pow(x2.val_, x1 - 1.0) / g;

      return fvar<T>(u, x2.d_ * der2);
    }
  }
}
#endif
