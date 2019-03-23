#ifndef STAN_MATH_REV_CORE_SCOPED_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_SCOPED_CHAINABLESTACK_HPP

#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/start_nested.hpp>
#include <stan/math/rev/core/recover_memory_nested.hpp>

namespace stan {
namespace math {

struct ScopedChainableStack {
  typedef ChainableStack::AutodiffStackStorage chainablestack_t;
  typedef ChainableStack::AutodiffStackQueue chainablequeue_t;
  typedef std::shared_ptr<chainablestack_t> stack_ptr_t;
  stack_ptr_t local_stack_;

  ScopedChainableStack()
      : local_stack_(ChainableStack::instance().get_child_stack()) {}

  ScopedChainableStack(chainablestack_t& parent_stack)
      : local_stack_(parent_stack.get_child_stack()) {}

  template <typename F>
  void execute(const F& f) {
    chainablequeue_t& local_queue = ChainableStack::queue();

    try {
      start_nested();
      const std::size_t nested_stack_instance = local_queue.current_instance_;
      stack_ptr_t nested_stack
          = local_queue.instance_stack_[nested_stack_instance];
      local_queue.instance_stack_[nested_stack_instance] = local_stack_;
      ChainableStack::instance_
          = local_queue.instance_stack_[nested_stack_instance].get();

      f();

      local_queue.instance_stack_[nested_stack_instance] = nested_stack;
      ChainableStack::instance_ = nested_stack.get();
      recover_memory_nested();
    } catch (const std::exception& e) {
      local_queue.instance_stack_[local_queue.current_instance_].reset(
          new chainablestack_t(ChainableStack::queue().stack_id_));
      recover_memory_nested();
      throw;
    }
  }

  void append_to_stack(chainablestack_t& destination_stack) {
    destination_stack.var_stack_.insert(destination_stack.var_stack_.end(),
                                        local_stack_->var_stack_.begin(),
                                        local_stack_->var_stack_.end());
    local_stack_->var_stack_.clear();
    destination_stack.var_nochain_stack_.insert(
        destination_stack.var_nochain_stack_.end(),
        local_stack_->var_nochain_stack_.begin(),
        local_stack_->var_nochain_stack_.end());
    local_stack_->var_nochain_stack_.clear();
    destination_stack.memalloc_.store_stack(local_stack_->memalloc_);
  }
};

}  // namespace math
}  // namespace stan
#endif
