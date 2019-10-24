#ifndef STAN_MATH_PRIM_SCAL_META_PROMOTE_ARGS_HPP
#define STAN_MATH_PRIM_SCAL_META_PROMOTE_ARGS_HPP

#include <boost/math/tools/promotion.hpp>

namespace stan {

/**
 * Convenience alias for boost tools promote_args
 */
template <typename... Args>
using promote_args_t = typename boost::math::tools::promote_args<Args...>::type;

}  // namespace stan
#endif
