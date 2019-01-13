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
  ChainableStack::instance().var_stack_.clear();
  ChainableStack::instance().var_nochain_stack_.clear();
  for (auto &x : ChainableStack::instance().var_alloc_stack_) {
    delete x;
  }
  ChainableStack::instance().var_alloc_stack_.clear();
  ChainableStack::instance().memalloc_.recover_all();
}

// recover memory the stack globally for all threads
static inline void recover_memory_global() {
  typedef ChainableStack::AutodiffStackStorage local_ad_stack_t;

  /*
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
  */
}

}  // namespace math
}  // namespace stan
#endif
