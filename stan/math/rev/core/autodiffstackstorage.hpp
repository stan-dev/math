#ifndef STAN_MATH_REV_CORE_AUTODIFFSTACKSTORAGE_HPP
#define STAN_MATH_REV_CORE_AUTODIFFSTACKSTORAGE_HPP

#include <stan/math/memory/stack_alloc.hpp>
#include <vector>
#include <tbb/concurrent_vector.h>

namespace stan {
namespace math {

/**
 * Provides a thread_local storage (TLS) singleton if needed. Read
 * warnings below!  For performance reasons the singleton is a global
 * static pointer which is initialized to a constant expression (the
 * nullptr). Doing allows for a faster access to the TLS pointer
 * instance. If one were to directly initialize the pointer with a
 * call to new, the compiler will insert additional code for each
 * reference in the code to the TLS. The additional code checks the
 * initialization status thus everytime one accesses the TLS, see [4]
 * for details and an example which demonstrates this. However, this
 * design requires that the pointer must be initialized with the
 * init() function before any var's are
 * instantiated. Otherwise a segfault will occur. The init() function
 * must be called for any thread which wants to use reverse mode
 * autodiff facilites.
 *
 * Furthermore, the declaration as a global (possibly thread local)
 * pointer allows the compiler to apply necessary inlining to get
 * maximal performance. However, this design suffers from "the static
 * init order fiasco"[0].  Anywhere this is used, we must be
 * absolutely positive that it doesn't matter when the singleton will
 * get initialized relative to other static variables.  In exchange,
 * we get a more performant singleton pattern.
 *
 * Formely, stan-math used in the threading case the defacto standard C++11
 * singleton pattern relying on a function wrapping a static local
 * variable. This standard pattern is expected to be well supported
 * by the major compilers (as its standard), but it does incur some
 * performance penalty.  There has been some discussion on this; see
 * [1] and [2] and the discussions those PRs link to as well. The
 * current design has a much reduced/almost no performance penalty
 * since the access to the TLS is not wrapped in extra function
 * calls. Moreover, the thread_local declaration below can be changed
 * to the __thread keyword which is a compiler-specific extension of
 * g++,clang++&intel and is around for much longer than C++11 such
 * that these compilers should support this design just as robust.
 *
 * These are thread_local only if the user asks for it with
 * -DSTAN_THREADS. This is primarily because Apple clang compilers
 * before 2016 don't support thread_local and the additional
 * performance cost. We have proposed removing support for those[3],
 * and at that time we should evaluate the performance of a switch to
 * thread_local.  If there is no loss in performance, we can remove
 * this ifdef.
 *
 * [0] https://isocpp.org/wiki/faq/ctors#static-init-order
 * [1] https://github.com/stan-dev/math/pull/840
 * [2] https://github.com/stan-dev/math/pull/826
 * [3]
 * http://discourse.mc-stan.org/t/potentially-dropping-support-for-older-versions-of-apples-version-of-clang/3780/
 * [4] https://discourse.mc-stan.org/t/thread-performance-penalty/7306/8
 */
template <typename ChainableT, typename ChainableAllocT>
struct AutodiffStackSingleton {
  typedef AutodiffStackSingleton<ChainableT, ChainableAllocT>
      AutodiffStackSingleton_t;

  struct AutodiffStackStorage {
    AutodiffStackStorage &operator=(const AutodiffStackStorage &) = delete;

    std::vector<ChainableT *> var_stack_;
    std::vector<ChainableT *> var_nochain_stack_;
    std::vector<ChainableAllocT *> var_alloc_stack_;
    stack_alloc memalloc_;

    // nested positions
    /*
    std::vector<size_t> nested_var_stack_sizes_;
    std::vector<size_t> nested_var_nochain_stack_sizes_;
    std::vector<size_t> nested_var_alloc_stack_starts_;
    */
  };

  struct AutodiffStackQueue {
    AutodiffStackQueue()
        : instance_stack_(1, std::shared_ptr<AutodiffStackStorage>(
                                 new AutodiffStackStorage())),
          current_instance_(0) {}

    AutodiffStackQueue &operator=(const AutodiffStackQueue &) = delete;

    std::vector<std::shared_ptr<AutodiffStackStorage>> instance_stack_;
    tbb::concurrent_vector<std::shared_ptr<AutodiffStackStorage>>
        instance_nochain_stack_;

    std::size_t current_instance_;
  };

  AutodiffStackSingleton() = delete;
  explicit AutodiffStackSingleton(AutodiffStackSingleton_t const &) = delete;
  AutodiffStackSingleton &operator=(const AutodiffStackSingleton_t &) = delete;

  constexpr static inline AutodiffStackStorage &instance() {
    return *instance_;
  }

  static inline AutodiffStackQueue &queue() {
    static
#ifdef STAN_THREADS
        thread_local
#endif
        AutodiffStackQueue queue_;
    return queue_;
  }

  static AutodiffStackStorage *init() {
    if (instance_ == nullptr) {
      AutodiffStackQueue &local_queue = queue();
      instance_
          = local_queue.instance_stack_[local_queue.current_instance_].get();
    }
    return instance_;
  }

  static
#ifdef STAN_THREADS
      thread_local
#endif
      AutodiffStackStorage *instance_;
};

template <typename ChainableT, typename ChainableAllocT>
#ifdef STAN_THREADS
thread_local
#endif
    typename AutodiffStackSingleton<ChainableT,
                                    ChainableAllocT>::AutodiffStackStorage
        *AutodiffStackSingleton<ChainableT, ChainableAllocT>::instance_
#ifdef STAN_THREADS
    = nullptr;
#else
    = AutodiffStackSingleton<ChainableT, ChainableAllocT>::init();
#endif

}  // namespace math
}  // namespace stan
#endif
