#ifndef STAN_MATH_REV_CORE_SCOPED_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_SCOPED_CHAINABLESTACK_HPP

#include <stan/math/rev/core/chainablestack.hpp>

#include <vector>

namespace stan {
namespace math {

/**
 * The AD tape of reverse mode AD is by default stored globally within the
 * process (or thread). With the ScopedChainableStack class one may execute a
 * nullary functor with reference to an AD tape which is stored with the
 * instance of ScopedChainableStack. Example:
 *
 * ScopedChainableStack scoped_stack;
 *
 * double cgrad = scoped_stack.execute([] {
 *   var a = 1;
 *   var b = 4;
 *   var c = a + b;
 *   grad();
 *   return c.adj();
 * });
 *
 * Doing so will not interfere with the process (or thread) AD tape.
 */
class ScopedChainableStack {
  ChainableStack::AutodiffStackStorage local_stack_;

  struct activate_scope {
    ChainableStack::AutodiffStackStorage& local_stack_;
    ChainableStack::AutodiffStackStorage* parent_stack_;

    explicit activate_scope(ChainableStack::AutodiffStackStorage& local_stack)
        : local_stack_(local_stack) {
      parent_stack_ = ChainableStack::instance_;
      ChainableStack::instance_ = &local_stack_;
    }

    ~activate_scope() {
      ChainableStack::instance_ = parent_stack_;
    }
  };

 public:
  ScopedChainableStack() = default;

  /**
   * ScopedChainableStack must call the destructor on all the
   * chainable_stack objects allocated on it
   */
  ~ScopedChainableStack() {
    for (size_t i = 0; i < local_stack_.var_alloc_stack_.size(); ++i) {
      delete local_stack_.var_alloc_stack_[i];
    }
  }

  /**
   * Execute in the current thread a nullary function and write the AD
   * tape to local_stack_ of this instance. The function may return
   * any type.
   *
   * @tparam F nullary functor to evaluate
   * @param f instance of nullary functor
   * @return Result of evaluated functor
   */
  template <typename F>
  auto execute(F&& f) {
    activate_scope active_scope(local_stack_);
    return std::forward<F>(f)();
  }
};

}  // namespace math
}  // namespace stan
#endif
