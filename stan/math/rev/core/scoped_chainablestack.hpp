#ifndef STAN_MATH_REV_CORE_SCOPED_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_SCOPED_CHAINABLESTACK_HPP

#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/start_nested.hpp>
#include <stan/math/rev/core/recover_memory_nested.hpp>

#include <vector>

namespace stan {
namespace math {

class ScopedChainableStack {
  ChainableStack::AutodiffStackStorage local_stack_;

  std::vector<ChainableStack::AutodiffStackStorage*> stack_queue_;

  struct activate_scope {
    ScopedChainableStack& scoped_stack_;

    explicit activate_scope(ScopedChainableStack& scoped_stack)
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
  ScopedChainableStack() = default;

  // execute in the current thread a nullary function and write the AD
  // tape to local_stack_ of this instance
  template <typename F>
  auto execute(F&& f) {
    activate_scope active_scope(*this);
    return std::forward<F>(f)();
  }
};

}  // namespace math
}  // namespace stan
#endif
