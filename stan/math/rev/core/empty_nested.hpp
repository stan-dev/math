#ifndef STAN_MATH_REV_CORE_EMPTY_NESTED_HPP
#define STAN_MATH_REV_CORE_EMPTY_NESTED_HPP

#include <stan/math/rev/core/chainablestack.hpp>

namespace stan {
namespace math {

/**
 * Return true if there is no nested autodiff being executed.
 */
static inline bool empty_nested() {
  // return ChainableStack::instance().nested_var_stack_sizes_.empty();
  // return ChainableStack::queue().current_instance_ == 0;
  // if the root is active, then we do not have nesting
  return ChainableStack::instance_->is_root();
}

}  // namespace math
}  // namespace stan
#endif
