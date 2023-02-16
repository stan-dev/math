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
 public:
  ScopedChainableStack() = default;

  /**
   * Execute in the current thread a function and write the AD
   * tape to local_stack_ of this instance. The function may return
   * any type.
   *
   * @tparam F functor to evaluate
   * @param f instance of functor
   * @param args arguments passed to functor
   * @return Result of evaluated functor
   */
  template <typename F, typename... Args>
  decltype(auto) execute(F&& f, Args&&... args) {
    nested_rev_autodiff nested;
    return std::forward<F>(f)(std::forward<Args>(args)...);
  }

  // Prevent undesirable operations
  ScopedChainableStack(const ScopedChainableStack&) = delete;
  ScopedChainableStack& operator=(const ScopedChainableStack&) = delete;
};

}  // namespace math
}  // namespace stan
#endif
