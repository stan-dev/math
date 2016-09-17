#ifndef STAN_MATH_PRIM_SCAL_FUN_BOOST_POLICY_HPP
#define STAN_MATH_PRIM_SCAL_FUN_BOOST_POLICY_HPP

#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>

namespace stan {
  namespace math {

    typedef boost::math::policies::policy<
      boost::math::policies::overflow_error<
        boost::math::policies::errno_on_error> >
    boost_policy_t;

  }
}
#endif
