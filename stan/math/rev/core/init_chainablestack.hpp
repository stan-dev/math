#ifndef STAN_MATH_REV_CORE_INIT_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_INIT_CHAINABLESTACK_HPP

#include <stan/math/rev/core/chainablestack.hpp>

namespace stan {
namespace math {
namespace {

/**
 * Initializes the AD stack for the main process. See
 * autodiffstackstorage.hpp for more explanations.
 *
 * TODO(wds15): remove once the Intel TBB is mandatory (and let
 * ad_tape_observer handle this).
 */
#ifndef STAN_THREADS
ChainableStack global_stack_instance_init;
#endif

}  // namespace
}  // namespace math
}  // namespace stan
#endif
