#ifndef STAN_MATH_PRIM_SCAL_FUN_ACOSH_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ACOSH_HPP

#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <boost/math/special_functions/acosh.hpp>

namespace stan {
  namespace math {

    /**
     * Return the inverse hyperbolic cosine of the specified value.
     *
     * @param[in] x Argument.
     * @return Inverse hyperbolic cosine of the argument.
     * @throw std::domain_error If argument is less than 1.
     */
    inline double acosh(double x) {
      return boost::math::acosh(x, boost_policy_t());
    }

    /**
     * Integer version of acosh.
     *
     * @param[in] x Argument.
     * @return Inverse hyperbolic cosine of the argument.
     * @throw std::domain_error If argument is less than 1.
     */
    inline double acosh(int x) {
      return acosh(static_cast<double>(x));
    }

  }
}
#endif
