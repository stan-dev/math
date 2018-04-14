#ifndef STAN_MATH_REV_CORE_AUTODIFFSTACKSTORAGE_HPP
#define STAN_MATH_REV_CORE_AUTODIFFSTACKSTORAGE_HPP

#include <stan/math/memory/stack_alloc.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename ChainableT, typename ChainableAllocT>
struct AutodiffStackStorage {
  std::vector<ChainableT*> var_stack_;
  std::vector<ChainableT*> var_nochain_stack_;
  std::vector<ChainableAllocT*> var_alloc_stack_;
  stack_alloc memalloc_;

  // nested positions
  std::vector<size_t> nested_var_stack_sizes_;
  std::vector<size_t> nested_var_nochain_stack_sizes_;
  std::vector<size_t> nested_var_alloc_stack_starts_;
};

/**
 * Provides a thread_local singleton. Read warnings below!
 * This is useful in a pretty limited set of circumstances, not just because
 * it's thread_local, but also because of "the static init order fiasco"[0].
 * Anywhere this is used, we must be absolutely positive that it doesn't matter
 * when the singleton will get initialized relative to other static variables.
 * In exchange, we get a much more performant singleton pattern than is
 * available with the standard C++11 singleton pattern relying on a function
 * wrapping a static local variable. There has been some discussion on this; see
 * [1] and [2] and the discussions those PRs link to as well.
 *
 * These are thread_local only if the user asks for it with -DSTAN_THREADS. This
 * is because Apple clang compilers before 2016 don't support thread_local. We
 * have proposed removing support for those[3], and will remove this ifdef at
 * that time.
 *
 * [0] https://isocpp.org/wiki/faq/ctors#static-init-order
 * [1] https://github.com/stan-dev/math/pull/840
 * [2] https://github.com/stan-dev/math/pull/826
 * [3]
 * http://discourse.mc-stan.org/t/potentially-dropping-support-for-older-versions-of-apples-version-of-clang/3780/
 */
template <typename ChainableT, typename ChainableAllocT>
struct AutodiffStackSingleton {
#ifdef STAN_THREADS
  thread_local
#endif
      static AutodiffStackStorage<ChainableT, ChainableAllocT>
          instance_;
};

template <typename ChainableT, typename ChainableAllocT>
#ifdef STAN_THREADS
thread_local
#endif
    AutodiffStackStorage<ChainableT, ChainableAllocT>
        AutodiffStackSingleton<ChainableT, ChainableAllocT>::instance_;

}  // namespace math
}  // namespace stan
#endif
