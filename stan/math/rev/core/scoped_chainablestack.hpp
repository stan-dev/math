#ifndef STAN_MATH_REV_CORE_SCOPED_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_SCOPED_CHAINABLESTACK_HPP

#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/start_nested.hpp>
#include <stan/math/rev/core/recover_memory_nested.hpp>

namespace stan {
namespace math {

class ScopedChainableStack {
  typedef ChainableStack::AutodiffStackStorage chainablestack_t;
  chainablestack_t& parent_stack_;
  chainablestack_t local_stack_;

 public:
  ScopedChainableStack() : parent_stack_(ChainableStack::instance()) {}

  ScopedChainableStack(chainablestack_t& parent_stack)
      : parent_stack_(parent_stack) {}

  // execute in the current thread a nullary function and write the AD
  // tape to local_stack_ of this instance
  template <typename F>
  void execute(const F& f) {
    if (local_stack_.is_active()) {
      f();
      return;
    }

    try {
      local_stack_.activate();
      f();
      local_stack_.deactivate();
    } catch (const std::exception& e) {
      local_stack_.deactivate();
      throw;
    }
  }

  void append_to_parent() {
    parent_stack_.var_stack_.insert(parent_stack_.var_stack_.end(),
                                    local_stack_.var_stack_.begin(),
                                    local_stack_.var_stack_.end());
    local_stack_.var_stack_.clear();
    parent_stack_.var_nochain_stack_.insert(
        parent_stack_.var_nochain_stack_.end(),
        local_stack_.var_nochain_stack_.begin(),
        local_stack_.var_nochain_stack_.end());
    local_stack_.var_nochain_stack_.clear();
    parent_stack_.memalloc_.store_stack(local_stack_.memalloc_);
  }

  void recover() { local_stack_.recover(); }
};

}  // namespace math
}  // namespace stan
#endif
