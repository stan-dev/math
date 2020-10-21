#ifndef STAN_MATH_REV_CORE_RECOVER_MEMORY_HPP
#define STAN_MATH_REV_CORE_RECOVER_MEMORY_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/empty_nested.hpp>
#include <stan/math/rev/core/needs_destructor.hpp>
#include <stdexcept>

namespace stan {
namespace math {

/**
 * Recover memory used for all variables for reuse.
 *
 * @throw std::logic_error if <code>empty_nested()</code> returns
 * <code>false</code>
 */
static inline void recover_memory() {
  if (!empty_nested()) {
    throw std::logic_error(
        "empty_nested() must be true"
        " before calling recover_memory()");
  }
  ChainableStack::instance_->var_stack_.clear();
  ChainableStack::instance_->var_nochain_stack_.clear();
  for (auto &x : ChainableStack::instance_->var_alloc_stack_) {
    delete x;
  }
  for (auto x : ChainableStack::instance_->destructor_stack_) {
    if(x!=nullptr){
      x->~needs_destructor();
    }
  }
  ChainableStack::instance_->var_alloc_stack_.clear();
  ChainableStack::instance_->memalloc_.recover_all();
}

}  // namespace math
}  // namespace stan
#endif
