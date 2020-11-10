#ifndef STAN_MATH_REV_CORE_SCOPED_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_SCOPED_CHAINABLESTACK_HPP

#include <stan/math/rev/core/chainablestack.hpp>

#include <mutex>
#include <stdexcept>
#include <thread>

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
 * double cgrad_a = scoped_stack.execute([] {
 *   var a = 2.0;
 *   var b = 4.0;
 *   var c = a*a + b;
 *   c.grad();
 *   return a.adj();
 * });
 *
 * Doing so will not interfere with the process (or thread) AD tape.
 */
class ScopedChainableStack {
  ChainableStack::AutodiffStackStorage local_stack_;
  std::mutex local_stack_mutex_;

  struct activate_scope {
    ChainableStack::AutodiffStackStorage* parent_stack_;
    ScopedChainableStack& scoped_stack_;

    explicit activate_scope(ScopedChainableStack& scoped_stack)
        : parent_stack_(ChainableStack::instance_),
          scoped_stack_(scoped_stack) {
      if (!scoped_stack_.local_stack_mutex_.try_lock()) {
        throw std::logic_error{"Cannot recurse same instance scoped stacks."};
      }
      ChainableStack::instance_ = &scoped_stack.local_stack_;
    }

    ~activate_scope() {
      scoped_stack_.local_stack_mutex_.unlock();
      ChainableStack::instance_ = parent_stack_;
    }
  };

 public:
  ScopedChainableStack() = default;

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
    activate_scope active_scope(*this);
    return std::forward<F>(f)();
  }
};

}  // namespace math
}  // namespace stan
#endif
