#ifndef STAN_MATH_REV_CORE_RECOVER_MEMORY_GLOBAL_HPP
#define STAN_MATH_REV_CORE_RECOVER_MEMORY_GLOBAL_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/empty_nested.hpp>
#include <stdexcept>

#include <iostream>

namespace stan {
namespace math {

/**
 * Recover memory global used for all variables for reuse in all
 * threads.
 *
 * @throw std::logic_error if <code>empty_nested()</code> returns
 * <code>false</code>
 */
static inline void recover_memory_global() {
  /* todo: check for each thread
  if (!empty_nested()) {
    throw std::logic_error(
        "empty_nested() must be true"
        " before calling recover_memory()");
  }
  */

  std::size_t i = 0;

  for (auto& kv : internal::global_observer.thread_tape_map_) {
    ChainableStack::AutodiffStackStorage* instance_ = kv.second->active_instance_;
    if (instance_ == nullptr) {
      std::cout << "nullptr for thread " << i << " ?? "<< std::endl;
      continue;
    }
    i++;
    instance_->var_stack_.clear();
    instance_->var_nochain_stack_.clear();
    for (auto& x : instance_->var_alloc_stack_) {
      delete x;
    }
    instance_->var_alloc_stack_.clear();
    instance_->memalloc_.recover_all();
  }
}

}  // namespace math
}  // namespace stan
#endif
