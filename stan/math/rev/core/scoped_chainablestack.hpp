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

  // execute in the current thread a nullary function and write the AD
  // tape to local_stack_ of this instance
  template <typename F>
  void execute(const F& f) {
    chainablequeue_t& local_queue = ChainableStack::queue();

    chainablestack_t* current_stack = ChainableStack::instance_;

    try {
      /*
      start_nested();
      // replace the nested AD tape with the local one from this
      // instance
      const std::size_t nested_stack_instance = local_queue.current_instance_;
      stack_ptr_t nested_stack
          = local_queue.instance_stack_[nested_stack_instance];
      local_queue.instance_stack_[nested_stack_instance] = local_stack_;
      ChainableStack::instance_
          = local_queue.instance_stack_[nested_stack_instance].get();

      f();

      // revert the replacement
      local_queue.instance_stack_[nested_stack_instance] = nested_stack;
      ChainableStack::instance_ = nested_stack.get();
      recover_memory_nested();
      */

      // slimmer version
      ChainableStack::instance_ = local_stack_.get();
      f();
      ChainableStack::instance_ = current_stack;

    } catch (const std::exception& e) {
      // local_queue.instance_stack_[local_queue.current_instance_].reset(
      //    new chainablestack_t(ChainableStack::queue().stack_id_));
      // recover_memory_nested();
      ChainableStack::instance_ = current_stack;
      throw;
    }
  }

  void append_to_stack(chainablestack_t& target_stack) {
    target_stack.var_stack_.insert(target_stack.var_stack_.end(),
                                   local_stack_->var_stack_.begin(),
                                   local_stack_->var_stack_.end());
    local_stack_->var_stack_.clear();
    target_stack.var_nochain_stack_.insert(
        target_stack.var_nochain_stack_.end(),
        local_stack_->var_nochain_stack_.begin(),
        local_stack_->var_nochain_stack_.end());
    local_stack_->var_nochain_stack_.clear();
    target_stack.memalloc_.store_stack(local_stack_->memalloc_);
  }

  void recover() {
    local_stack_->var_stack_.clear();
    local_stack_->var_nochain_stack_.clear();

    for (size_t i = 0; i < local_stack_->var_alloc_stack_.size(); ++i) {
      delete local_stack_->var_alloc_stack_[i];
    }
    local_stack_->var_alloc_stack_.clear();
    local_stack_->memalloc_.recover_all();
  }
};

}  // namespace math
}  // namespace stan
#endif
