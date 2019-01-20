#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>

namespace stan {
namespace math {

class vari;
class chainable_alloc;

typedef AutodiffStackSingleton<vari, chainable_alloc> ChainableStack;

/**
 * Instantiates an instance of the ChainableStack if not already
 * initialized. This function must be called before any autodiff
 * variable get's instantiated within any thread which performs
 * reverse mode autodiff operations.
 */
static inline ChainableStack::AutodiffStackStorage* init() {
  if (ChainableStack::instance_ == nullptr) {
    ChainableStack::AutodiffStackQueue& queue = ChainableStack::queue();
    ChainableStack::instance_ = queue.instance_stack_[0].get();
  }
  return ChainableStack::instance_;
}

}  // namespace math
}  // namespace stan
#endif
