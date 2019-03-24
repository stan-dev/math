#ifndef STAN_MATH_REV_CORE_INIT_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_INIT_CHAINABLESTACK_HPP

#include <stan/math/rev/core/chainablestack.hpp>

namespace stan {
namespace math {
namespace {

/**
 * Initializes the AD stack for the main process. See
 * autodiffstackstorage.hpp for more explanations.
 */
ChainableStack global_stack_instance_init;
}  // namespace
}  // namespace math
}  // namespace stan
#endif
