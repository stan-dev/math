#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_INST_CPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_INST_CPP

#include <stan/math/rev/core/chainablestack.hpp>

namespace stan {
namespace math {
namespace internal {

const ChainableStack::AutodiffStackStorage* __chainable_stack
    = ChainableStack::instantiate();
}
}  // namespace math
}  // namespace stan
#endif
