#ifndef STAN_MATH_REV_CORE_RECOVER_MEMORY_HPP
#define STAN_MATH_REV_CORE_RECOVER_MEMORY_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/empty_nested.hpp>
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
  if (!empty_nested())
    throw std::logic_error(
        "empty_nested() must be true"
        " before calling recover_memory()");

  /*
  auto clean_instance = [](ChainableStack::AutodiffStackStorage& instance) {
    instance.var_stack_.clear();
    instance.var_nochain_stack_.clear();
    for (auto& x : instance.var_alloc_stack_) {
      delete x;
    }
    instance.var_alloc_stack_.clear();
    instance.memalloc_.recover_all();
  };

  clean_instance(ChainableStack::instance());

  for (auto& instance_stack_ptr : ChainableStack::queue().instance_stack_) {
    clean_instance(*instance_stack_ptr);
  }
  */

  ChainableStack::instance().recover();

  // not needed, since this function may only be called if we are not nested
  // ChainableStack::queue().current_instance_ = 0;

  /*
  const std::size_t stack_id = ChainableStack::queue().stack_id_;

  std::lock_guard<std::mutex> global_stack_lock(
      ChainableStack::global_stack_mutex_);

  ChainableStack::global_stack_t old_global_stack;
  ChainableStack::global_stack_t& global_stack = ChainableStack::global_stack();

  old_global_stack.swap(global_stack);

  global_stack.reserve(old_global_stack.size());

  for (auto& instance_ptr : old_global_stack) {
    if (instance_ptr) {
      if (instance_ptr->stack_id_ == stack_id) {
        clean_instance(*instance_ptr);
        instance_ptr.reset();
      } else {
        global_stack.emplace_back(instance_ptr);
      }
    }
  }
  */
}

// recover memory the stack globally for all threads

/*
static inline void recover_memory_global() {
typedef ChainableStack::AutodiffStackStorage local_ad_stack_t;

std::for_each(ChainableStack::instance_.begin(),
              ChainableStack::instance_.end(),
              [](local_ad_stack_t &local_instance) {
                if (!local_instance.nested_var_stack_sizes_.empty())
                  throw std::logic_error(
                      "empty_nested() must be true"
                      " before calling recover_memory_global()");
                local_instance.var_stack_.clear();
                local_instance.var_nochain_stack_.clear();
                for (auto &x : local_instance.var_alloc_stack_) {
                  delete x;
                }
                local_instance.var_alloc_stack_.clear();
                local_instance.memalloc_.recover_all();
              });
}
*/

}  // namespace math
}  // namespace stan
#endif
