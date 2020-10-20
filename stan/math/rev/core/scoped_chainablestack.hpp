#ifndef STAN_MATH_REV_CORE_SCOPED_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_SCOPED_CHAINABLESTACK_HPP

#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/start_nested.hpp>
#include <stan/math/rev/core/recover_memory_nested.hpp>

#include <vector>

namespace stan {
namespace math {

class ScopedChainableStack {
  typedef ChainableStack::AutodiffStackStorage chainablestack_t;
  chainablestack_t local_stack_;

  using stack_queue_t = std::vector<chainablestack_t*>;
  stack_queue_t stack_queue_;

  struct activate_scope {
    ScopedChainableStack& scoped_stack_;

    activate_scope(ScopedChainableStack& scoped_stack)
        : scoped_stack_(scoped_stack) {
      scoped_stack_.stack_queue_.push_back(ChainableStack::instance_);
      ChainableStack::instance_ = &scoped_stack_.local_stack_;
    }

    ~activate_scope() {
      ChainableStack::instance_ = scoped_stack_.stack_queue_.back();
      scoped_stack_.stack_queue_.pop_back();
    }
  };

 public:
  ScopedChainableStack() {}

  // execute in the current thread a nullary function and write the AD
  // tape to local_stack_ of this instance
  template <typename F>
  F execute(F f) {
    activate_scope active_scope(*this);
    f();

    return f;
  }

  /*
  template <typename F>
  void execute(const F& f) {
    // It's actually impossible to leave the stack in an active state
    // behind, but if that happens, then we may not try-catch to
    // ensure we deactivate.
    if (local_stack_.is_active()) {
      f();
      //return f;
    }

    try {
      local_stack_.activate();
      f();
      local_stack_.deactivate();
    } catch (const std::exception& e) {
      local_stack_.deactivate();
      throw;
    }
    //return f;
  }
  */
};

}  // namespace math
}  // namespace stan
#endif
