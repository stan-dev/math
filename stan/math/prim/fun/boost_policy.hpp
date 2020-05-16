#ifndef STAN_MATH_PRIM_FUN_BOOST_POLICY_HPP
#define STAN_MATH_PRIM_FUN_BOOST_POLICY_HPP

#include <stan/math/prim/meta.hpp>
#include <boost/math/policies/policy.hpp>

namespace stan {
namespace math {

/**
 * Boost policy that overrides the defaults to match the built-in
 * C++ standard library functions.
 *
 * The non-default behavior from Boost's built-ins are
 * (1) overflow errors return error numbers on error.
 * (2) pole errors return error numbers on error.
 * (3) doubles passed to Boost are not promoted to long double.
 *
 * The policy is equipped with an optional generic argument B controlling the
 * precision in some functions. If set to 0, the maximum precision available
 * in the type being used is demanded from Boost. Otherwise, it correspond to
 * the approximately B-bit precision, i.e. for trading speed for accuracy.
 */
template <int B = 0>
using boost_policy_t = boost::math::policies::policy<
    boost::math::policies::overflow_error<
        boost::math::policies::errno_on_error>,
    boost::math::policies::pole_error<boost::math::policies::errno_on_error>,
    boost::math::policies::promote_double<false>,
    boost::math::policies::digits2<B>>;
}  // namespace math
}  // namespace stan
#endif
